"""
Gestionnaire centralisé des données et calculs
"""

import numpy as np
import pyvista as pv
import glob
import os
import heapq
from scipy.spatial import cKDTree, ConvexHull
from scipy.stats import skew as scipy_skew, kurtosis as scipy_kurtosis
from typing import Dict, Optional, Tuple


class DataManager:
    """Charge et prépare toutes les données du patient"""
    
    def __init__(self, folder_path: str):
        self.path = os.path.normpath(folder_path)
        self.meshes = {}
        self.computed_metrics = {}
        
        # Paramètres globaux
        self.SEARCH_RADIUS = 4.0
        self.MIN_T = 1.0
        
        self.CV_HEALTHY = 0.6   # v_max (Cedilnik EP-Europace 2018)
        self.CV_BORDER = 0.3    # v(p) — vitesse au point d'inflexion
        self.CV_SCAR = 0.01      # référence informative (non utilisé dans la simu sigmoïde)
        
        # Flag pour marquer si cicatrice estimée (vs IRM)
        self._has_estimated_scar = False
        
    # ------------------------------------------------------------------
    # Méthodes utilitaires internes
    # ------------------------------------------------------------------

    @staticmethod
    def _clip_annular_rings(mesh: pv.PolyData, threshold: float = 0.5,
                            dilation: int = 2) -> pv.PolyData:
        """Retire les anneaux mitral/aortique (dist EPI < threshold) du mesh EPI DIST MAP.

        La zone retirée est élargie de `dilation` rangées de voisins afin de
        gratter un peu plus de rayon autour de l'anneau, puis on ne conserve que
        la plus grande composante connexe — ce qui élimine au passage les îlots
        résiduels parfois laissés au centre de l'anneau."""
        if mesh.active_scalars is None or mesh.n_points == 0:
            return mesh

        scalar_name = mesh.active_scalars_name
        scalars = mesh.active_scalars.astype(float)

        # Identifier les points annulaires (distance EPI ≈ 0)
        near_zero_pts = set(int(i) for i in np.where(scalars < threshold)[0])
        if not near_zero_pts:
            return mesh

        faces = mesh.regular_faces  # (n_cells, 3)

        # Dilatation : tout point appartenant à une face qui touche déjà la zone
        # annulaire est ajouté à la sélection. En répétant, on élargit le rayon
        # retiré de quelques rangées de triangles autour de l'anneau.
        for _ in range(max(0, dilation)):
            sel = np.fromiter(near_zero_pts, dtype=int)
            touching = (
                np.isin(faces[:, 0], sel) |
                np.isin(faces[:, 1], sel) |
                np.isin(faces[:, 2], sel)
            )
            grown = np.unique(faces[touching].ravel())
            if len(grown) == len(near_zero_pts):
                break
            near_zero_pts.update(int(i) for i in grown)

        sel = np.fromiter(near_zero_pts, dtype=int)
        keep_mask = ~(
            np.isin(faces[:, 0], sel) |
            np.isin(faces[:, 1], sel) |
            np.isin(faces[:, 2], sel)
        )

        clean_indices = np.where(keep_mask)[0]
        if len(clean_indices) == 0:
            return mesh

        mesh_clean = mesh.extract_cells(clean_indices)
        mesh_poly = mesh_clean.extract_surface()

        mesh_conn = mesh_poly.connectivity(extraction_mode='largest')
        # Connectivity écrase le scalaire actif avec RegionId → le restaurer
        if scalar_name and scalar_name in mesh_conn.point_data.keys():
            mesh_conn.set_active_scalars(scalar_name)

        return mesh_conn

    # ------------------------------------------------------------------

    def load_patient_data(self) -> bool:
        """Charge tous les fichiers VTK disponibles"""
        
        file_patterns = {
            'lv_epi_dist': '*LV EPI DIST MAP*.vtk',
            # 'lv_epi_dist': '*LV ENDO DIST MAP*.vtk',
            'lv_endo_dist': '*LV ENDO DIST MAP*.vtk',
            'dense_scar': '*DENSE SCAR*.vtk',
            'scar': '*SCAR (LE)*.vtk',
            'scar_endo': '*SCAR_ENDO*.vtk',
            'scar_epi': '*SCAR_EPI*.vtk',
            'scar_intra': '*SCAR_INTRA*.vtk',
            'scar_transmurality': '*SCAR TRANSMURALITY*.vtk',
            'lv_fat': '*LV FAT*.vtk',
        }
        
        for key, pattern in file_patterns.items():
            files = glob.glob(os.path.join(self.path, "**", pattern), recursive=True)
            if files:
                self.meshes[key] = pv.read(files[0])
        
        # Retirer automatiquement les anneaux mitral/aortique du EPI DIST MAP
        if 'lv_epi_dist' in self.meshes:
            self.meshes['lv_epi_dist'] = self._clip_annular_rings(self.meshes['lv_epi_dist'])
        
        # Pré-calculer les métriques de base si les données sont là
        if 'lv_epi_dist' in self.meshes:
            self.compute_wall_thickness()
            self.compute_ciaccio_ratio()
            self.compute_deceleration_zones()
            self.compute_channels()
            self.compute_isthmus_prediction()
            self.compute_laplacian()
            self.compute_scar_burden()
            self.compute_channelness()
            self.compute_cv_map()
            self.compute_border_zone_cedilnik()
            self.compute_combined_zones()
        if 'scar_transmurality' in self.meshes:
            self.compute_transmurality()
        
        # Si IRM manquant, estimer la cicatrice depuis l'épaisseur
        if 'dense_scar' not in self.meshes and 'wall_thickness' in self.computed_metrics:
            self.estimate_scar_from_thickness()
        
        return len(self.meshes) > 0
    
    def get_dense_scar(self) -> Optional[pv.PolyData]:
        """Retourne le mesh de cicatrice dense"""
        return self.meshes.get('dense_scar')
    
    def get_scar_distribution(self) -> Optional[pv.PolyData]:
        """Retourne le mesh de distribution de cicatrice par localisation"""
        # On retourne scar_endo, scar_epi, scar_intra combinés
        parts = []
        for key in ['scar_endo', 'scar_epi', 'scar_intra']:
            if key in self.meshes:
                parts.append(self.meshes[key])
        if not parts:
            return None
        return parts[0] if len(parts) == 1 else parts[0].merge(parts[1:])
    
    def get_fat_mesh(self) -> Optional[pv.PolyData]:
        """Retourne le mesh de graisse"""
        return self.meshes.get('lv_fat')
    
    def compute_wall_thickness(self) -> Optional[pv.PolyData]:
        """Calcul de l'épaisseur de paroi.
        
        La carte LV EPI DIST MAP mesure, en chaque point de la surface EPI,
        la distance jusqu'à la surface ENDO la plus proche — c'est directement
        l'épaisseur pariétale. On l'utilise telle quelle, sans ajout de la distance ENDO.
        """
        if 'lv_epi_dist' not in self.meshes:
            return None
        
        mesh_epi = self.meshes['lv_epi_dist']
        epi = mesh_epi.active_scalars.astype(float)
        
        mesh = mesh_epi.copy()
        mesh["Wall_Thickness"] = epi   # EPI dist MAP = épaisseur pariétale réelle
        mesh["EPI_Distance"]   = epi   # alias conservé pour compatibilité
        self.computed_metrics['wall_thickness'] = mesh

        # Auto-détecter les seuils cliniques (2/5 mm) et stocker les stats patient
        self.detect_scar_thresholds()

        return mesh

    def detect_scar_thresholds(self) -> dict:
        """Détecte automatiquement les seuils de cicatrice d'après l'épaisseur de paroi.

        Valeurs cliniques de référence :
          - Dense scar : T < 2 mm  (Stevenson 1989 — conduction block, voltage < 0.5 mV)
          - Scar total : T < 5 mm  (Cedilnik EP-Europace 2018 — zone altérée)
          - Tissu sain : T ≥ 5 mm

        Ces seuils calibrent automatiquement le point d'inflexion sigmoïde p
        pour compute_channelness : p = (2+5)/2 = 3.5 mm.
        """
        scar_mm       = 5.0    # Cedilnik 2018
        dense_scar_mm = 2.0    # Stevenson 1989
        self._scar_thresh       = scar_mm
        self._dense_scar_thresh = dense_scar_mm

        result = {
            'scar_thresh_mm':       scar_mm,
            'dense_scar_thresh_mm': dense_scar_mm,
            'p_inflection_mm':      (scar_mm + dense_scar_mm) / 2.0,  # 3.5 mm
            'method': 'clinical_standard',
            'ref': 'Cedilnik EP-Europace 2018; Stevenson 1989',
        }

        if 'wall_thickness' in self.computed_metrics:
            T = self.computed_metrics['wall_thickness']['Wall_Thickness'].astype(float)
            n = len(T)
            result['n_pts']          = n
            result['wt_mean_mm']     = float(np.mean(T))
            result['wt_median_mm']   = float(np.median(T))
            result['wt_p5_mm']       = float(np.percentile(T, 5))
            result['pct_dense_scar'] = float(np.mean(T < dense_scar_mm) * 100)
            result['pct_scar']       = float(np.mean((T >= dense_scar_mm) & (T < scar_mm)) * 100)
            result['pct_healthy']    = float(np.mean(T >= scar_mm) * 100)

        self._scar_threshold_stats = result
        return result
    
    def compute_ciaccio_ratio(self) -> Optional[pv.PolyData]:
        """Calcul du ratio de courbure Ciaccio"""
        
        if 'wall_thickness' not in self.computed_metrics:
            self.compute_wall_thickness()
        
        mesh = self.computed_metrics['wall_thickness'].copy()
        pts = mesh.points
        T = mesh["Wall_Thickness"].astype(float)
        
        tree = cKDTree(pts)
        neighbors = tree.query_ball_point(pts, r=self.SEARCH_RADIUS)
        
        rho = np.zeros(len(T))
        
        for i, idx in enumerate(neighbors):
            if len(idx) < 2 or T[i] < self.MIN_T:
                continue
            Tn = T[idx]
            delta_T = np.max(Tn) - T[i]
            rho[i] = delta_T / (self.SEARCH_RADIUS * T[i])
        
        mesh["Ciaccio_Ratio"] = rho
        mesh["Ciaccio_Ratio_Display"] = np.clip(rho, 0, 0.5)
        
        self.computed_metrics['ciaccio_ratio'] = mesh
        
        return mesh
    
    def compute_tri(self) -> Optional[pv.PolyData]:
        """Terrain Ruggedness Index (TRI).
        
        Pour chaque point, calcule la racine carrée de la somme des carrés
        des différences d'épaisseur avec ses voisins directs (adjacents par arêtes).
        
        TRI = sqrt(sum((T_neighbor - T_central)^2))
        
        Un TRI élevé indique une surface très irrégulière (rugueuse),
        caractéristique des zones de cicatrice hétérogène ou des canaux.
        """
        if 'wall_thickness' not in self.computed_metrics:
            self.compute_wall_thickness()
        
        mesh = self.computed_metrics['wall_thickness'].copy()
        T = mesh["Wall_Thickness"].astype(float)
        n = len(T)
        
        # Construire le graphe d'adjacence
        edges = mesh.extract_all_edges()
        adjacency = [[] for _ in range(n)]
        for i in range(edges.n_cells):
            cell = edges.get_cell(i)
            if cell.n_points == 2:
                p1, p2 = cell.point_ids
                adjacency[p1].append(p2)
                adjacency[p2].append(p1)
        
        # Calculer TRI pour chaque point
        tri = np.zeros(n)
        for i in range(n):
            neighbors = adjacency[i]
            if len(neighbors) == 0:
                continue
            
            # Différences avec les voisins
            diffs = T[neighbors] - T[i]
            
            # TRI = sqrt(sum(diff^2))
            tri[i] = np.sqrt(np.sum(diffs ** 2))
        
        mesh["TRI"] = tri
        tri_p95 = np.percentile(tri, 95) if np.any(tri > 0) else 1.0
        mesh["TRI_Normalized"] = np.clip(tri / max(tri_p95, 1e-6), 0, 1)
        self.computed_metrics['tri'] = mesh
        return mesh
    
    def compute_deceleration_zones(self, dz_threshold: float = 0.33) -> Optional[pv.PolyData]:
        """Zones de décélération avec seuil paramétrable"""
        
        if 'ciaccio_ratio' not in self.computed_metrics:
            self.compute_ciaccio_ratio()
        
        mesh = self.computed_metrics['ciaccio_ratio'].copy()
        rho = mesh["Ciaccio_Ratio"]
        
        mesh["DZ_Mask"] = (rho >= dz_threshold).astype(float)
        self.computed_metrics['deceleration_zones'] = mesh
        
        return mesh
    
    def compute_channels(self, target_thickness: float = None) -> Optional[pv.PolyData]:
        """Canaux de conduction (0-5mm).

        Optionnellement calcule un champ `Channel_Target_Display` centré
        sur `target_thickness` (mm) pour mise en évidence par un gradient
        rouge symétrique autour de l'épaisseur cible.
        """
        
        if 'wall_thickness' not in self.computed_metrics:
            self.compute_wall_thickness()
        
        mesh = self.computed_metrics['wall_thickness'].copy()
        T = mesh["Wall_Thickness"].astype(float)
        pts = mesh.points

        # Canaux : épaisseur entre 0 et 5mm
        channel_mask = (T <= 5.0).astype(bool)

        # Region simplifiée : 0 = reste, 1 = canal
        region = np.zeros(len(T), dtype=float)
        region[channel_mask] = 1.0

        mesh["Channel_Zone"] = channel_mask.astype(float)
        mesh["Channel_Region"] = region
        mesh["Channel_Score"] = np.where(channel_mask, 1 - (np.abs(T - 2.5) / 2.5), 0)

        # --- Target thickness display (highlight centered on target_thickness) ---
        # Si `target_thickness` fourni, construire un champ [0,1] où 1 = épaisseur cible
        # et décroît symétriquement vers 0. On utilise une tolérance douce pour
        # créer un dégradé (tolérance par défaut 0.5 mm).
        if target_thickness is not None:
            tol = 0.5  # mm, demi-largeur de la bande de mise en évidence
            diff = np.abs(T - float(target_thickness))
            intensity = np.clip(1.0 - diff / tol, 0.0, 1.0)
            mesh["Channel_Target_Display"] = intensity.astype(float)
            try:
                mesh.field_data['channel_target'] = np.array([float(target_thickness)])
            except Exception:
                pass

        self.computed_metrics['channels'] = mesh
        return mesh

    def update_channel_target(self, target_thickness: float, tol: float = 0.5) -> Optional[pv.PolyData]:
        """Met à jour en place le champ `Channel_Target_Display` pour les channels existants.

        Retourne le mesh mis à jour ou None si channels non calculés.
        """
        if 'channels' not in self.computed_metrics:
            return None

        mesh = self.computed_metrics['channels']
        if 'Wall_Thickness' not in mesh.array_names:
            return None

        T = mesh['Wall_Thickness'].astype(float)
        diff = np.abs(T - float(target_thickness))
        intensity = np.clip(1.0 - diff / float(tol), 0.0, 1.0)

        # Remplacer ou créer le tableau en place
        try:
            mesh.point_data['Channel_Target_Display'] = intensity.astype(float)
        except Exception:
            mesh['Channel_Target_Display'] = intensity.astype(float)

        try:
            mesh.field_data['channel_target'] = np.array([float(target_thickness)])
        except Exception:
            pass

        return mesh
        
        self.computed_metrics['channels'] = mesh
        
        return mesh
    
    def compute_laplacian(self) -> Optional[pv.PolyData]:
        """Calcul du Laplacien de l'épaisseur de paroi (divergence du gradient).
        
        Le Laplacien met en évidence :
        - Valeurs positives : zones convexes / creux (amincissement local)
        - Valeurs négatives : zones concaves / bosses (épaississement local)
        - ~ 0 : zones planes / uniformes
        
        Utilise PyVista compute_derivative pour le gradient puis divergence.
        """
        if 'wall_thickness' not in self.computed_metrics:
            self.compute_wall_thickness()
        
        mesh = self.computed_metrics['wall_thickness'].copy()
        T = mesh["Wall_Thickness"].astype(float)
        
        # Assigner le scalaire actif pour le calcul du gradient
        mesh["_T_for_grad"] = T
        
        try:
            # 1. Calculer le gradient de l'épaisseur (vecteur 3D)
            grad_mesh = mesh.compute_derivative(
                scalars="_T_for_grad",
                gradient=True,
                faster=False
            )
            
            # 2. Calculer la divergence du gradient = Laplacien
            # Le gradient est un champ vectoriel (3 composantes)
            gradient = grad_mesh["gradient"]
            
            # Ajouter les composantes du gradient séparément pour calculer la divergence
            grad_mesh["grad_x"] = gradient[:, 0]
            grad_mesh["grad_y"] = gradient[:, 1]
            grad_mesh["grad_z"] = gradient[:, 2]
            
            # Divergence = d(grad_x)/dx + d(grad_y)/dy + d(grad_z)/dz
            div_x = grad_mesh.compute_derivative(scalars="grad_x", gradient=True, faster=False)
            div_y = grad_mesh.compute_derivative(scalars="grad_y", gradient=True, faster=False)
            div_z = grad_mesh.compute_derivative(scalars="grad_z", gradient=True, faster=False)
            
            # Extraire les composantes diagonales du Hessien
            laplacian = (div_x["gradient"][:, 0] + 
                        div_y["gradient"][:, 1] + 
                        div_z["gradient"][:, 2])
            
            # Clipper les valeurs extrêmes (artefacts aux bords du mesh)
            p1, p99 = np.percentile(laplacian, [1, 99])
            laplacian = np.clip(laplacian, p1, p99)
            
            mesh["Laplacian"] = laplacian
            mesh["Laplacian_Abs"] = np.abs(laplacian)
            
            self.computed_metrics['laplacian'] = mesh
            return mesh
            
        except Exception as e:
            print(f"Laplacian computation error: {e}")
            import traceback
            traceback.print_exc()
            return None
    
    def compute_isthmus_prediction(self) -> Optional[pv.PolyData]:
        """Prédiction isthmes VT depuis topographie CT (Takigawa, Heart Rhythm 2019).
        Score combiné : canal 1–5 mm + proximité scar + gradient + étroitesse."""
        if 'wall_thickness' not in self.computed_metrics:
            self.compute_wall_thickness()
        if 'wall_thickness' not in self.computed_metrics:
            return None
        
        mesh = self.computed_metrics['wall_thickness'].copy()
        T = mesh["Wall_Thickness"].astype(float)
        pts = mesh.points
        n = len(T)
        
        # ── Critère 1 : Zone candidate (1–5 mm, Takigawa : >1mm) ──
        candidate_mask = (T >= 1.0) & (T < 5.0)
        
        # ── Critère 2 : Proximité du scar dense (<2 mm) ──
        dense_scar_ids = np.where(T < 2.0)[0]
        proximity_to_scar = np.zeros(n)
        if len(dense_scar_ids) > 0:
            tree_scar = cKDTree(pts[dense_scar_ids])
            dist_to_scar, _ = tree_scar.query(pts, k=1)
            proximity_to_scar = np.exp(-dist_to_scar / 5.0)
        
        # ── Critère 3 : Gradient d'épaisseur (adjacence directe) ──
        edges = mesh.extract_all_edges()
        adjacency = [[] for _ in range(n)]
        for i in range(edges.n_cells):
            cell = edges.get_cell(i)
            if cell.n_points == 2:
                p1, p2 = cell.point_ids
                adjacency[p1].append(p2)
                adjacency[p2].append(p1)
        
        gradient_score = np.zeros(n)
        for i in range(n):
            if not candidate_mask[i] or len(adjacency[i]) == 0:
                continue
            neighbor_T = T[adjacency[i]]
            gradient_score[i] = max(T[i] - np.min(neighbor_T), 0) / max(T[i], 0.1)
        
        # ── Critère 4 : Étroitesse locale du canal ──
        # Distance aux points non-candidats → diamètre approximatif
        non_cand_ids = np.where(~candidate_mask | (T < 1.0))[0]
        channel_width = np.full(n, 50.0)
        if len(non_cand_ids) > 0:
            tree_nc = cKDTree(pts[non_cand_ids])
            dist_nc, _ = tree_nc.query(pts, k=1)
            channel_width = dist_nc * 2.0
        
        narrowness_score = np.clip(1.0 - channel_width / 30.0, 0, 1)
        
        # ── Score combiné ──
        isthmus_score = np.zeros(n)
        isthmus_score[candidate_mask] = (
            0.35 * proximity_to_scar[candidate_mask] +
            0.30 * gradient_score[candidate_mask] +
            0.35 * narrowness_score[candidate_mask]
        )
        
        max_score = np.max(isthmus_score) if np.max(isthmus_score) > 0 else 1.0
        isthmus_score = isthmus_score / max_score
        
        # Zones critiques (Takigawa : centre WT ~2.4mm, bords 1.5mm)
        isthmus_critical = (isthmus_score > 0.5) & (T >= 1.5) & (T <= 4.0)
        
        mesh["Isthmus_Probability"] = isthmus_score
        mesh["Isthmus_Critical"] = isthmus_critical.astype(float)
        mesh["Isthmus_Proximity_Scar"] = proximity_to_scar
        mesh["Isthmus_Width_mm"] = channel_width
        mesh["Isthmus_Boundaries"] = (gradient_score > 0.3).astype(float)
        
        self._isthmus_stats = {
            'n_isthmus_pts': int(np.sum(isthmus_score > 0.3)),
            'n_critical_pts': int(np.sum(isthmus_critical)),
            'isthmus_pct': float(np.mean(isthmus_score > 0.3) * 100),
            'max_score': float(np.max(isthmus_score)),
            'mean_width_mm': float(np.mean(channel_width[candidate_mask])) if np.any(candidate_mask) else 0.0,
        }
        self.computed_metrics['isthmus'] = mesh
        return mesh
    
    def compute_scar_burden(self, thresh_healthy=5.0, thresh_border=4.0, thresh_dense=2.0) -> Optional[pv.PolyData]:
        """
        Classification de l'épaisseur de paroi en 4 zones (seuils configurables) :
          0 = Healthy   : T > thresh_healthy
          1 = Border    : thresh_border < T ≤ thresh_healthy
          2 = Scar      : thresh_dense < T ≤ thresh_border
          3 = Dense     : T ≤ thresh_dense

        Utah grade calculé sur le % surface (scar+dense = T≤thresh_healthy) :
          I < 5%  /  II 5-20%  /  III 20-35%  /  IV > 35%
        """
        if 'wall_thickness' not in self.computed_metrics:
            self.compute_wall_thickness()
        if 'wall_thickness' not in self.computed_metrics:
            return None

        mesh = self.computed_metrics['wall_thickness'].copy()
        T = mesh["Wall_Thickness"].astype(float)

        # --- Classification 4 zones (seuils configurables) ---
        scar_zone = np.zeros(len(T), dtype=int)
        scar_zone[T <= thresh_healthy] = 1
        scar_zone[T <= thresh_border] = 2
        scar_zone[T <= thresh_dense] = 3
        scar_zone[T > thresh_healthy] = 0

        mesh["Scar_Zone"] = scar_zone.astype(float)
        mesh["Scar_Burden_Display"] = scar_zone.astype(float)

        # Store thresholds in field_data for the renderer
        mesh.field_data['thresh_healthy'] = np.array([thresh_healthy])
        mesh.field_data['thresh_border'] = np.array([thresh_border])
        mesh.field_data['thresh_dense'] = np.array([thresh_dense])

        # --- Surface pondérée par cellule ---
        cell_mesh = mesh.point_data_to_cell_data()
        cell_sizes = cell_mesh.compute_cell_sizes()
        cells_T    = cell_sizes.cell_data['Wall_Thickness']
        cell_areas = cell_sizes.cell_data['Area']
        total_area = float(np.sum(cell_areas))

        healthy_area = float(np.sum(cell_areas[cells_T > thresh_healthy]))
        border_area  = float(np.sum(cell_areas[(cells_T > thresh_border) & (cells_T <= thresh_healthy)]))
        scar_area    = float(np.sum(cell_areas[(cells_T > thresh_dense) & (cells_T <= thresh_border)]))
        dense_area   = float(np.sum(cell_areas[cells_T <= thresh_dense]))

        def _pct(a): return (a / total_area * 100) if total_area > 0 else 0.0
        healthy_pct = _pct(healthy_area)
        border_pct  = _pct(border_area)
        scar_pct    = _pct(scar_area)
        dense_pct   = _pct(dense_area)
        # Utah calculé sur T ≤ 5mm (border+scar+dense)
        utah_pct = border_pct + scar_pct + dense_pct

        if utah_pct < 5.0:
            utah_grade = 'I'
            utah_label = 'I - Minimal'
        elif utah_pct < 20.0:
            utah_grade = 'II'
            utah_label = 'II - Mild'
        elif utah_pct < 35.0:
            utah_grade = 'III'
            utah_label = 'III - Moderate'
        else:
            utah_grade = 'IV'
            utah_label = 'IV - Severe'

        # field_data pour le renderer et le rapport
        mesh.field_data['wt_healthy_pct'] = np.array([healthy_pct])
        mesh.field_data['wt_border_pct']  = np.array([border_pct])
        mesh.field_data['wt_scar_pct']    = np.array([scar_pct])
        mesh.field_data['wt_dense_pct']   = np.array([dense_pct])
        mesh.field_data['utah_label']     = np.array([utah_label], dtype=object)
        mesh.field_data['total_area_cm2'] = np.array([total_area / 100.0])
        # rétro-compatibilité
        mesh.field_data['scar_burden_pct']  = np.array([utah_pct])
        mesh.field_data['dense_burden_pct'] = np.array([dense_pct])
        mesh.field_data['utah_grade']       = np.array([ord(utah_grade[0])])

        self._scar_burden = {
            'utah_label':    utah_label,
            'utah_pct':      utah_pct,
            'healthy_pct':   healthy_pct,
            'border_pct':    border_pct,
            'scar_pct':      scar_pct,
            'dense_pct':     dense_pct,
            'total_area_cm2': total_area / 100.0,
        }

        self.computed_metrics['scar_burden'] = mesh
        return mesh
    
    def compute_transmurality(self) -> Optional[pv.PolyData]:
        """Analyse de la transmuralité"""
        
        if 'scar_transmurality' not in self.meshes:
            return None
        
        mesh = self.meshes['scar_transmurality'].copy()
        
        if mesh.active_scalars is not None:
            trans = mesh.active_scalars
        elif len(mesh.point_data) > 0:
            trans = mesh.point_data[list(mesh.point_data.keys())[0]]
        else:
            return None
        
        mesh["Transmurality"] = trans
        self.computed_metrics['transmurality'] = mesh
        
        return mesh
    
    def compute_channelness(self, n_pacing=30, v_max=0.6, v_min=0.01, p=3.0, r=2.0) -> Optional[pv.PolyData]:
        """Channelness map via simulation Eikonal multi-site (Cedilnik, EP-Europace 2018).
        Fonction de transfert sigmoïde T→CV, Dijkstra sur graphe du maillage.
        Paramètres : n_pacing (sites), v_max/v_min (m/s), p (inflexion mm), r (raideur)."""
        if 'wall_thickness' not in self.computed_metrics:
            self.compute_wall_thickness()
        if 'wall_thickness' not in self.computed_metrics:
            return None
        
        mesh = self.computed_metrics['wall_thickness'].copy()
        T = mesh["Wall_Thickness"].astype(float)
        pts = mesh.points
        n = len(T)
        
        # Step 1: Cedilnik logistic transfer function (T → CV)
        cv = (v_max - v_min) / (1.0 + np.exp(r * (p - T))) + v_min
        mesh["CV_Cedilnik"] = cv
        
        # Tissue classification
        tissue_class = np.zeros(n, dtype=float)
        tissue_class[T < 2.0] = 2.0                        # Dense scar
        tissue_class[(T >= 2.0) & (T < 5.0)] = 1.0         # Gray zone / channel candidates
        mesh["Tissue_Class"] = tissue_class
        
        # Step 2: Build mesh adjacency graph
        edges = mesh.extract_all_edges()
        adjacency = [[] for _ in range(n)]
        for i in range(edges.n_cells):
            cell = edges.get_cell(i)
            if cell.n_points == 2:
                p1, p2 = cell.point_ids
                d = float(np.linalg.norm(pts[p1] - pts[p2]))
                adjacency[p1].append((p2, d))
                adjacency[p2].append((p1, d))
        
        # Step 3: Select pacing sites (farthest-point sampling)
        healthy_ids = np.where(T >= 5.0)[0]
        if len(healthy_ids) < max(n_pacing, 2):
            healthy_ids = np.argsort(T)[-max(n_pacing * 2, 20):]
        
        # Start with the thickest point, then greedily pick farthest
        selected = [int(healthy_ids[np.argmax(T[healthy_ids])])]
        for _ in range(min(n_pacing, len(healthy_ids)) - 1):
            dists = np.min(
                np.linalg.norm(
                    pts[healthy_ids, None, :] - pts[selected][None, :, :],
                    axis=2),
                axis=1)
            selected.append(int(healthy_ids[np.argmax(dists)]))
        pacing_sites = selected
        
        # Step 4: Dijkstra simulations (pathological vs healthy reference)
        def _dijkstra_mesh(start_idx, velocities):
            """Dijkstra shortest-time on mesh with edge-averaged velocities."""
            act = np.full(n, np.inf)
            act[start_idx] = 0.0
            heap = [(0.0, int(start_idx))]
            visited = np.zeros(n, dtype=bool)
            while heap:
                t_cur, u = heapq.heappop(heap)
                if visited[u]:
                    continue
                visited[u] = True
                for v_idx, d in adjacency[u]:
                    if visited[v_idx]:
                        continue
                    avg_v = (velocities[u] + velocities[v_idx]) * 0.5
                    if avg_v < 1e-4:
                        continue
                    new_t = t_cur + d / avg_v
                    if new_t < act[v_idx]:
                        act[v_idx] = new_t
                        heapq.heappush(heap, (new_t, int(v_idx)))
            return act
        
        cv_ref = np.full(n, v_max)  # Uniform healthy reference
        
        delay_sum = np.zeros(n, dtype=float)
        delay_count = np.zeros(n, dtype=int)
        
        for start in pacing_sites:
            act_patho = _dijkstra_mesh(start, cv)
            act_ref   = _dijkstra_mesh(start, cv_ref)
            
            valid = (~np.isinf(act_patho)) & (~np.isinf(act_ref)) & (act_ref > 1e-6)
            ratio = np.ones(n, dtype=float)
            ratio[valid] = act_patho[valid] / act_ref[valid]
            delay_sum[valid] += ratio[valid]
            delay_count[valid] += 1
        
        delay_count = np.maximum(delay_count, 1)
        mean_delay = delay_sum / delay_count
        mesh["Activation_Delay"] = mean_delay
        
        # Step 5: Full-ventricle channelness score
        # Channels (2–5 mm, high delay) → brightest; dense scar → dimmer; healthy → dark.
        
        mask_channel = (T >= 2.0) & (T < 5.0)
        mask_scar    = T < 2.0
        mask_healthy = T >= 5.0
        
        channelness = np.zeros(n, dtype=float)
        
        # (a) Healthy tissue → near zero (no isthmus)
        if np.any(mask_healthy):
            channelness[mask_healthy] = np.clip(
                0.05 * (1.0 - (T[mask_healthy] - 5.0) / 10.0), 0.0, 0.05)
        
        # (b) Dense scar → moderate (visible substrate, not a channel)
        if np.any(mask_scar):
            scar_severity = np.clip(1.0 - T[mask_scar] / 2.0, 0.0, 1.0)
            channelness[mask_scar] = 0.15 + 0.15 * scar_severity  # 0.15–0.30
        
        # (c) Channel zone (2–5 mm) → scored by simulation delay
        if np.any(mask_channel):
            ch_delay = mean_delay[mask_channel]
            # Normalize within channel band (robust scaling)
            ch_p95 = np.percentile(ch_delay, 95) if len(ch_delay) > 0 else 1.0
            ch_base = np.min(ch_delay) if len(ch_delay) > 0 else 1.0
            denom = max(ch_p95 - ch_base, 1e-6)
            ch_norm = np.clip((ch_delay - ch_base) / denom, 0.0, 1.0)
            channelness[mask_channel] = 0.40 + 0.60 * ch_norm  # 0.40–1.00
        
        # (d) Smooth transitions (5–7 mm)
        mask_trans = (T >= 5.0) & (T <= 7.0)
        if np.any(mask_trans):
            alpha = np.clip((7.0 - T[mask_trans]) / 2.0, 0.0, 1.0)
            trans_delay = mean_delay[mask_trans]
            trans_base = np.min(trans_delay) if len(trans_delay) > 0 else 1.0
            trans_norm = np.clip((trans_delay - trans_base) / max(denom, 1e-6), 0.0, 0.5)
            channelness[mask_trans] = alpha * (0.20 + 0.20 * trans_norm)
        
        mesh["Channelness"] = channelness
        
        # Step 6: Statistics
        channel_threshold = 0.5
        n_channel_pts = int(np.sum(channelness >= channel_threshold))
        max_ch = float(np.max(channelness))
        ch_active = channelness[channelness > 0.1]
        mean_ch = float(np.mean(ch_active)) if len(ch_active) > 0 else 0.0
        
        mean_delay_ch = float(np.mean(mean_delay[mask_channel])) if np.any(mask_channel) else 1.0
        max_delay_ch  = float(np.max(mean_delay[mask_channel]))  if np.any(mask_channel) else 1.0
        
        self._channelness_stats = {
            'n_channel_points': n_channel_pts,
            'n_channel_pct': float(n_channel_pts / n * 100) if n > 0 else 0.0,
            'max_channelness': max_ch,
            'mean_channelness': mean_ch,
            'mean_delay_ratio': mean_delay_ch,
            'max_delay_ratio': max_delay_ch,
            'n_pacing_sites': len(pacing_sites),
            'cv_range_ms': (float(np.min(cv)), float(np.max(cv))),
            'transfer_params': {'v_max': v_max, 'v_min': v_min, 'p': p, 'r': r},
        }
        
        mesh.field_data['channelness_max'] = np.array([max_ch])
        mesh.field_data['mean_delay'] = np.array([mean_delay_ch])
        mesh.field_data['n_pacing'] = np.array([len(pacing_sites)])
        mesh.field_data['p_threshold'] = np.array([p])
        
        self.computed_metrics['channelness'] = mesh
        return mesh

    def compute_anatomical_channelness(
        self,
        h_min: float = 3.0,
        max_width: float = 8.0,
    ) -> Optional[pv.PolyData]:
        """Geometric channelness from wall thickness alone — no electrophysiology.

        Algorithm
        ---------
        1. Viable mask  : T > h_min
        2. Distance transform : d(x) = dist to nearest non-viable point (cKDTree)
           → local corridor width w(x) = 2·d(x)
        3. Narrow corridor mask : viable AND w(x) < max_width
        4. Skeleton (medial axis) : viable points that are local maxima of d(x)
           among mesh-edge neighbours, restricted to the narrow corridor
        5. Connected components of narrow corridors (union-find on mesh adjacency)
        6. Per-component L / W ratio
           L = 2 × max(‖skel_pt − centroid‖) along skeleton points
           W = mean corridor width of the component
        7. Combined score = 0.6 × C_point(x) + 0.4 × LW_norm(x)
           where C_point = clip(1 − w(x)/max_width, 0, 1)

        Args:
            h_min     : viable threshold (mm). Points with T ≤ h_min are non-viable.
            max_width : max corridor width (mm). Wider regions are healthy, not channels.
        """
        from collections import defaultdict

        if 'lv_epi_dist' not in self.meshes:
            return None

        mesh = self.meshes['lv_epi_dist'].copy()

        # Ensure thickness is present
        if 'EPI_Distance' not in mesh.array_names:
            thick = self.compute_wall_thickness()
            if thick is None:
                return None
            if 'EPI_Distance' not in thick.array_names:
                return None
            mesh = thick

        T   = mesh['EPI_Distance'].astype(float)
        pts = mesh.points
        n   = len(T)

        # ── Step 1 : viable / non-viable ─────────────────────────────────
        viable     = T > h_min
        non_viable = ~viable

        viable_ids     = np.where(viable)[0]
        non_viable_ids = np.where(non_viable)[0]

        if len(viable_ids) == 0:
            return None

        # ── Step 2 : distance transform → local corridor width ───────────
        anat_half_width = np.zeros(n)
        if len(non_viable_ids) > 0:
            nv_tree          = cKDTree(pts[non_viable_ids])
            d_vals, _        = nv_tree.query(pts[viable_ids], k=1)
            anat_half_width[viable_ids] = d_vals
        else:
            # Entire mesh is viable: use T as proxy
            anat_half_width = T / 2.0

        anat_width = anat_half_width * 2.0          # diameter
        mesh['Anat_Width'] = anat_width

        # ── Step 3 : narrow corridor mask ────────────────────────────────
        narrow_mask = viable & (anat_width > 0) & (anat_width < max_width)
        mesh['Narrow_Mask'] = narrow_mask.astype(float)

        # ── Step 4 : mesh adjacency (viable points only) + skeleton ──────
        adjacency = [[] for _ in range(n)]
        try:
            edges = mesh.extract_all_edges()
            for k in range(edges.n_cells):
                cell = edges.get_cell(k)
                if cell.n_points == 2:
                    p1, p2 = cell.point_ids
                    if viable[p1] and viable[p2]:
                        adjacency[p1].append(p2)
                        adjacency[p2].append(p1)
        except Exception:
            pass

        is_skeleton = np.zeros(n, dtype=bool)
        for vid in viable_ids:
            nbrs = adjacency[vid]
            if not nbrs:
                is_skeleton[vid] = True
                continue
            if anat_half_width[vid] >= max(anat_half_width[nb] for nb in nbrs):
                is_skeleton[vid] = True

        is_skeleton = is_skeleton & narrow_mask
        mesh['Is_Skeleton'] = is_skeleton.astype(float)

        # ── Step 5 : connected components (union-find on narrow mask) ─────
        narrow_ids  = np.where(narrow_mask)[0]
        narrow_set  = set(narrow_ids.tolist())
        g2l         = {gid: lid for lid, gid in enumerate(narrow_ids)}

        parent = list(range(len(narrow_ids)))

        def _find(x: int) -> int:
            root = x
            while parent[root] != root:
                root = parent[root]
            while parent[x] != root:
                parent[x], x = root, parent[x]
            return root

        for gid in narrow_ids:
            lid = g2l[gid]
            for nb in adjacency[gid]:
                if nb in narrow_set:
                    r1, r2 = _find(lid), _find(g2l[nb])
                    if r1 != r2:
                        parent[r1] = r2

        components: dict = defaultdict(list)
        for lid in range(len(narrow_ids)):
            components[_find(lid)].append(narrow_ids[lid])

        # ── Step 6 : L / W ratio per component ──────────────────────────
        lw_score   = np.zeros(n)
        lw_raw_max = 0.0

        for comp_gids in components.values():
            if len(comp_gids) < 4:
                continue
            comp_arr  = np.array(comp_gids)
            comp_pts  = pts[comp_arr]

            # Skeleton subset of this component
            skel_mask = is_skeleton[comp_arr]
            skel_pts  = comp_pts[skel_mask]

            if len(skel_pts) >= 2:
                centroid = skel_pts.mean(axis=0)
                L = 2.0 * float(np.max(np.linalg.norm(skel_pts - centroid, axis=1)))
            else:
                bb_min = comp_pts.min(axis=0)
                bb_max = comp_pts.max(axis=0)
                L = float(np.linalg.norm(bb_max - bb_min))

            W_mean = float(np.mean(anat_width[comp_arr]))
            if W_mean < 1e-3:
                continue

            lw = L / W_mean
            lw_score[comp_arr] = lw
            lw_raw_max = max(lw_raw_max, lw)

        # Normalise L/W → [0, 1] using p95
        lw_valid = lw_score[lw_score > 0]
        if len(lw_valid) > 0:
            p95 = float(np.percentile(lw_valid, 95)) if len(lw_valid) >= 20 \
                  else float(lw_valid.max())
            lw_score = np.clip(lw_score / max(p95, 1.0), 0.0, 1.0)

        mesh['LW_Ratio'] = lw_score

        # ── Step 7 : pointwise inverse-width score ────────────────────────
        c_point = np.zeros(n)
        c_point[viable_ids] = np.clip(
            1.0 - anat_width[viable_ids] / max(max_width, 1.0), 0.0, 1.0)

        # ── Step 8 : combined anatomical channelness ──────────────────────
        anat_ch = np.zeros(n)
        anat_ch[viable_ids] = (0.6 * c_point[viable_ids]
                               + 0.4 * lw_score[viable_ids])
        anat_ch = np.clip(anat_ch, 0.0, 1.0)
        mesh['Anatomical_Channelness'] = anat_ch

        # ── Field data for renderer annotations ──────────────────────────
        narrow_ch  = anat_ch[narrow_mask]
        viable_pct = float(np.mean(viable) * 100)
        narrow_pct = float(np.mean(narrow_mask) * 100)

        mesh.field_data['h_min']          = np.array([h_min])
        mesh.field_data['max_width']      = np.array([max_width])
        mesh.field_data['viable_pct']     = np.array([viable_pct])
        mesh.field_data['narrow_pct']     = np.array([narrow_pct])
        mesh.field_data['n_skeleton_pts'] = np.array([int(np.sum(is_skeleton))])
        mesh.field_data['mean_anat_ch']   = np.array(
            [float(np.mean(narrow_ch))] if len(narrow_ch) > 0 else [0.0])
        mesh.field_data['max_lw_raw']     = np.array([lw_raw_max])

        self.computed_metrics['anatomical_channelness'] = mesh
        return mesh

    def compute_cv_map(self, v_max=0.6, v_min=0.01, p=3.0, r=2.0) -> Optional[pv.PolyData]:
        """Carte CV via la fonction de transfert sigmoïde de Cedilnik (EP-Europace 2018).
        v(x) = (v_max − v_min) / (1 + exp(r·(p − W(x)))) + v_min"""
        if 'wall_thickness' not in self.computed_metrics:
            self.compute_wall_thickness()
        if 'wall_thickness' not in self.computed_metrics:
            return None
        
        mesh = self.computed_metrics['wall_thickness'].copy()
        T = mesh["Wall_Thickness"].astype(float)
        
        # Eq. 3 (thèse §2.2.2.5)
        cv = (v_max - v_min) / (1.0 + np.exp(r * (p - T))) + v_min
        
        mesh["CV_ms"] = cv
        mesh["CV_Normalized"] = cv / v_max
        
        # Classification fonctionnelle
        func_class = np.zeros(len(T), dtype=float)
        func_class[cv >= 0.3 * v_max] = 2.0    # Normal
        func_class[(cv >= 0.05) & (cv < 0.3 * v_max)] = 1.0  # Lente
        # 0.0 = bloquée
        mesh["CV_Functional_Class"] = func_class
        
        self._cv_stats = {
            'cv_mean': float(np.mean(cv)),
            'cv_std': float(np.std(cv)),
            'cv_min': float(np.min(cv)),
            'cv_max': float(np.max(cv)),
            'pct_blocked': float(np.mean(cv < 0.05) * 100),
            'pct_slow': float(np.mean((cv >= 0.05) & (cv < 0.3 * v_max)) * 100),
            'pct_normal': float(np.mean(cv >= 0.3 * v_max) * 100),
            'transfer_params': {'v_max': v_max, 'v_min': v_min, 'p': p, 'r': r},
        }
        
        self.computed_metrics['cv_map'] = mesh
        return mesh
    
    def compute_border_zone_cedilnik(self, dilation_radius: float = 2.0) -> Optional[pv.PolyData]:
        """Zone bordante Cedilnik : intersection de la dilatation scar ET tissu sain
        (dilation_radius=2mm par défaut). Retourne Cedilnik_BZ, Cedilnik_BZ_Mask."""
        if 'wall_thickness' not in self.computed_metrics:
            self.compute_wall_thickness()
        if 'wall_thickness' not in self.computed_metrics:
            return None
        
        mesh = self.computed_metrics['wall_thickness'].copy()
        T = mesh["Wall_Thickness"].astype(float)
        pts = mesh.points
        n = len(T)
        
        scar_mask = T < 5.0
        healthy_mask = T >= 5.0
        
        scar_ids = np.where(scar_mask)[0]
        healthy_ids = np.where(healthy_mask)[0]
        
        if len(scar_ids) == 0 or len(healthy_ids) == 0:
            # Pas de frontière à définir
            mesh["Cedilnik_BZ"] = np.zeros(n, dtype=float)
            mesh["Cedilnik_Distance"] = np.zeros(n, dtype=float)
            self.computed_metrics['border_zone_cedilnik'] = mesh
            return mesh
        
        # --- Dilatation morphologique ---
        tree_all = cKDTree(pts)
        
        # 1. Dilated scar : points à ≤ radius d'au moins un point scar
        tree_scar = cKDTree(pts[scar_ids])
        dist_to_scar, _ = tree_scar.query(pts, k=1)
        dilated_scar = dist_to_scar <= dilation_radius
        
        # 2. Dilated healthy : points à ≤ radius d'au moins un point sain
        tree_healthy = cKDTree(pts[healthy_ids])
        dist_to_healthy, _ = tree_healthy.query(pts, k=1)
        dilated_healthy = dist_to_healthy <= dilation_radius
        
        # 3. Border zone = intersection
        border_zone = dilated_scar & dilated_healthy
        
        # --- Classification ---
        # 0 = hors zone (ni scar ni healthy proche)
        # 1 = border zone (Cedilnik)
        # 2 = scar pur (dans scar, pas dans dilated_healthy)
        # 3 = sain pur (dans healthy, pas dans dilated_scar)
        classification = np.zeros(n, dtype=int)
        classification[scar_mask & ~border_zone] = 2   # scar pur
        classification[healthy_mask & ~border_zone] = 3 # sain pur
        classification[border_zone] = 1                  # border zone
        
        # Distance au bord scar/healthy le plus proche (pour gradient de visualisation)
        cedilnik_dist = np.minimum(dist_to_scar, dist_to_healthy)
        cedilnik_dist[~border_zone] = 0.0
        
        mesh["Cedilnik_BZ"] = classification.astype(float)
        mesh["Cedilnik_Distance"] = cedilnik_dist
        mesh["Cedilnik_BZ_Mask"] = border_zone.astype(float)
        
        # Stats
        bz_count = int(np.sum(border_zone))
        self._cedilnik_stats = {
            'bz_point_count': bz_count,
            'bz_point_pct': float(bz_count / n * 100),
            'dilation_radius_mm': dilation_radius,
        }
        
        # Surface de la BZ pondérée par cellule
        try:
            cell_mesh = mesh.point_data_to_cell_data()
            cell_sizes = cell_mesh.compute_cell_sizes()
            cells_bz = cell_sizes.cell_data['Cedilnik_BZ_Mask']
            cell_areas = cell_sizes.cell_data['Area']
            bz_area = float(np.sum(cell_areas[cells_bz > 0.5]))
            total_area = float(np.sum(cell_areas))
            self._cedilnik_stats['bz_area_cm2'] = bz_area / 100.0
            self._cedilnik_stats['bz_area_pct'] = (bz_area / total_area * 100) if total_area > 0 else 0.0
        except:
            pass
        
        self.computed_metrics['border_zone_cedilnik'] = mesh
        return mesh
    
    def compute_combined_zones(self) -> Optional[pv.PolyData]:
        """Zones combinées multi-critères : canal étroit + gradient + BZ Cedilnik + channelness."""
        # S'assurer que toutes les métriques sont calculées
        if 'wall_thickness' not in self.computed_metrics:
            self.compute_wall_thickness()
        if 'channelness' not in self.computed_metrics:
            self.compute_channelness()
        if 'border_zone_cedilnik' not in self.computed_metrics:
            self.compute_border_zone_cedilnik()
        if 'ciaccio_ratio' not in self.computed_metrics:
            self.compute_ciaccio_ratio()
        
        mesh = self.computed_metrics['wall_thickness'].copy()
        T = mesh["Wall_Thickness"].astype(float)
        n = len(T)
        
        # --- Critère 1 : Canal étroit (1-4mm) ---
        # Plus critique que 1-5mm, centré sur les vrais canaux de conduction
        narrow_channel = ((T >= 1.0) & (T <= 4.0)).astype(float)
        
        # --- Critère 2 : Gradient d'épaisseur élevé ---
        gradient_score = np.zeros(n, dtype=float)
        if 'ciaccio_ratio' in self.computed_metrics:
            rho = self.computed_metrics['ciaccio_ratio']["Ciaccio_Ratio"]
            # Normaliser [0, 1] avec saturation à 0.5
            gradient_score = np.clip(rho / 0.5, 0, 1)
        
        # --- Critère 3 : Border zone Cedilnik ---
        bz_score = np.zeros(n, dtype=float)
        if 'border_zone_cedilnik' in self.computed_metrics:
            bz_mesh = self.computed_metrics['border_zone_cedilnik']
            bz_score = bz_mesh["Cedilnik_BZ_Mask"].astype(float)
        
        # --- Critère 4 : Channelness ---
        ch_score = np.zeros(n, dtype=float)
        if 'channelness' in self.computed_metrics:
            ch_score = self.computed_metrics['channelness']["Channelness"].astype(float)
        
        # --- Score composite ---
        # Pondérations basées sur la pertinence arythmogène
        W_CHANNEL = 0.35    # Canal étroit : haute importance
        W_GRADIENT = 0.20   # Gradient : modéré
        W_BZ = 0.20         # Border zone Cedilnik  
        W_CHANNELNESS = 0.25 # Mesure morphologique
        
        combined = (W_CHANNEL * narrow_channel + 
                   W_GRADIENT * gradient_score +
                   W_BZ * bz_score +
                   W_CHANNELNESS * ch_score)
        
        # Normaliser [0, 1]
        if np.max(combined) > 0:
            combined = combined / np.max(combined)
        
        # --- Classification discrète ---
        # 0 = safe, 1 = low risk, 2 = moderate, 3 = high risk, 4 = critical
        zone = np.zeros(n, dtype=int)
        zone[combined >= 0.2] = 1
        zone[combined >= 0.4] = 2
        zone[combined >= 0.6] = 3
        zone[combined >= 0.8] = 4
        
        mesh["Combined_Score"] = combined
        mesh["Combined_Zone"] = zone.astype(float)
        
        self.computed_metrics['combined_zones'] = mesh
        return mesh

    # Métriques combinables disponibles pour le score custom
    COMBINABLE_METRICS = {
        'narrow_channel': {
            'label': 'Canal étroit (WT 1-4mm)',
            'description': 'Zone WT entre 1 et 4mm — couloir de conduction critique',
        },
        'wall_thickness': {
            'label': 'Épaisseur de paroi (classique)',
            'description': 'Carte d\'épaisseur pariétale — paroi fine = substrat',
        },
        'laplacian': {
            'label': 'Laplacien (épaisseur)',
            'description': 'Laplacien de l\'épaisseur — fortes variations locales de paroi',
        },
        'ciaccio': {
            'label': 'Gradient (Ciaccio)',
            'description': 'Changement brutal d\'épaisseur — front d\'onde ralenti',
        },
        'channelness': {
            'label': 'Channelness (Cedilnik)',
            'description': 'Simulation Eikonal — probabilité de canal de conduction',
        },
        'isthmus': {
            'label': 'Isthme (Takigawa)',
            'description': 'Probabilité d\'isthme de TV',
        },
        'deceleration': {
            'label': 'Zones de décélération',
            'description': 'Zones où le front d\'onde ralentit (ρ > seuil)',
        },
        'border_zone': {
            'label': 'Border zone (Cedilnik)',
            'description': 'Zone bordante dilatée (intersection scar+healthy) — Cedilnik 2018',
        },
        'local_entropy': {
            'label': 'Entropie locale',
            'description': 'Irrégularité morphologique locale — hétérogénéité du substrat',
        },
        'tri': {
            'label': 'TRI (Rugosité terrain)',
            'description': 'Terrain Ruggedness Index — variation d\'épaisseur avec voisins',
        },
        'cv_slow': {
            'label': 'Conduction lente (CV inversée)',
            'description': 'Vitesse de conduction faible — zones à risque de bloc',
        },
        'scar_proximity': {
            'label': 'Proximité du scar dense',
            'description': 'Score de proximité géographique au scar dense (< 2mm) — Stevenson 1989',
        },
        'activation_rv_late': {
            'label': 'Activation tardive — Pacing VD',
            'description': 'Zones à activation retardée lors d\'un pacing VD (normalisé 0→1)',
        },
        'activation_lv_late': {
            'label': 'Activation tardive — Pacing VG',
            'description': 'Zones à activation retardée lors d\'un pacing VG (normalisé 0→1)',
        },
    }

    def compute_custom_combined(self, config: Optional[list] = None) -> Optional[pv.PolyData]:
        """Score composite configurable — somme pondérée de métriques normalisées [0,1].
        config : [{'metric': str, 'weight': float}, ...] (clés dans COMBINABLE_METRICS).
        Retourne le mesh avec Combined_Custom [0,1]."""
        if 'wall_thickness' not in self.computed_metrics:
            self.compute_wall_thickness()
        if 'wall_thickness' not in self.computed_metrics:
            return None

        # Config par défaut (isthme critique = combinaison 1 de la littérature)
        if config is None or len(config) == 0:
            config = [
                {'metric': 'narrow_channel', 'weight': 0.35},
                {'metric': 'channelness',    'weight': 0.30},
                {'metric': 'ciaccio',        'weight': 0.20},
                {'metric': 'border_zone',    'weight': 0.15},
            ]

        mesh = self.computed_metrics['wall_thickness'].copy()
        T    = mesh["Wall_Thickness"].astype(float)
        n    = len(T)

        # S'assurer que les métriques dépendantes sont calculées
        if 'ciaccio_ratio' not in self.computed_metrics:
            self.compute_ciaccio_ratio()
        if 'channelness' not in self.computed_metrics:
            self.compute_channelness()
        if 'isthmus' not in self.computed_metrics:
            self.compute_isthmus_prediction()
        if 'deceleration_zones' not in self.computed_metrics:
            self.compute_deceleration_zones()
        if 'border_zone_cedilnik' not in self.computed_metrics:
            self.compute_border_zone_cedilnik()
        if 'local_entropy' not in self.computed_metrics:
            try:
                self.compute_local_entropy()
            except Exception:
                pass
        if 'tri' not in self.computed_metrics:
            try:
                self.compute_tri()
            except Exception:
                pass
        if 'cv_map' not in self.computed_metrics:
            try:
                self.compute_cv_map()
            except Exception:
                pass
        if 'laplacian' not in self.computed_metrics:
            try:
                self.compute_laplacian()
            except Exception:
                pass

        def _get_layer(metric_key: str) -> np.ndarray:
            """Extrait et normalise une couche [0, 1]."""
            arr = np.zeros(n, dtype=float)
            try:
                if metric_key == 'narrow_channel':
                    arr = ((T >= 1.0) & (T <= 4.0)).astype(float)
                elif metric_key == 'wall_thickness':
                    # Carte d'épaisseur classique : paroi fine = substrat → score
                    # élevé pour les parois amincies (inverse de l'épaisseur).
                    p99 = np.percentile(T, 99) if T.max() > 0 else 1.0
                    if p99 > 0:
                        arr = np.clip(1.0 - T / p99, 0, 1)
                elif metric_key == 'laplacian':
                    if 'laplacian' in self.computed_metrics:
                        lap = np.abs(self.computed_metrics['laplacian']["Laplacian"].astype(float))
                        p99 = np.percentile(lap, 99) if lap.max() > 0 else 1.0
                        arr = np.clip(lap / p99, 0, 1) if p99 > 0 else lap
                elif metric_key == 'ciaccio':
                    if 'ciaccio_ratio' in self.computed_metrics:
                        rho = self.computed_metrics['ciaccio_ratio']["Ciaccio_Ratio"].astype(float)
                        arr = np.clip(rho / 0.5, 0, 1)
                elif metric_key == 'channelness':
                    if 'channelness' in self.computed_metrics:
                        arr = np.clip(self.computed_metrics['channelness']["Channelness"].astype(float), 0, 1)
                elif metric_key == 'isthmus':
                    if 'isthmus' in self.computed_metrics:
                        arr = np.clip(self.computed_metrics['isthmus']["Isthmus_Probability"].astype(float), 0, 1)
                elif metric_key == 'deceleration':
                    if 'deceleration_zones' in self.computed_metrics:
                        arr = self.computed_metrics['deceleration_zones']["DZ_Mask"].astype(float)
                elif metric_key == 'border_zone':
                    if 'border_zone_cedilnik' in self.computed_metrics:
                        arr = self.computed_metrics['border_zone_cedilnik']["Cedilnik_BZ_Mask"].astype(float)
                elif metric_key == 'local_entropy':
                    if 'local_entropy' in self.computed_metrics:
                        v = self.computed_metrics['local_entropy']["Local_Entropy"].astype(float)
                        mx = v.max()
                        arr = v / mx if mx > 0 else v
                elif metric_key == 'tri':
                    if 'tri' in self.computed_metrics:
                        v = self.computed_metrics['tri']["TRI"].astype(float)
                        p99 = np.percentile(v, 99) if v.max() > 0 else 1.0
                        arr = np.clip(v / p99, 0, 1)
                elif metric_key == 'cv_slow':
                    if 'cv_map' in self.computed_metrics:
                        cv = self.computed_metrics['cv_map']["CV_ms"].astype(float)
                        mx = cv.max()
                        arr = 1.0 - (cv / mx) if mx > 0 else np.zeros(n)
                elif metric_key == 'scar_proximity':
                    # Distance au point le plus proche avec T < 2mm
                    dense_pts = mesh.points[T <= 2.0]
                    if len(dense_pts) > 0:
                        tree_dense = cKDTree(dense_pts)
                        dists, _ = tree_dense.query(mesh.points, k=1)
                        # Score élevé à proximité (< 10mm), décroît avec distance
                        arr = np.clip(1.0 - dists / 10.0, 0, 1)
                elif metric_key in ('activation_rv_late', 'activation_lv_late'):
                    # Activation tardive = zones qui s'activent en dernier
                    # → score élevé = candidats à ralentissement / re-entry
                    pacing = 'RV' if metric_key == 'activation_rv_late' else 'LV'
                    cache = f'activation_{pacing}'
                    if cache not in self.computed_metrics:
                        self.simulate_activation(pacing_site=pacing)
                    if cache in self.computed_metrics:
                        act = self.computed_metrics[cache]["Activation_Time"].astype(float)
                        # Interpoler sur le mesh courant si tailles différentes
                        src = self.computed_metrics[cache]
                        if len(act) == n:
                            arr_act = act
                        else:
                            tree_src = cKDTree(src.points)
                            _, idx = tree_src.query(mesh.points, k=1)
                            arr_act = act[idx]
                        mx = arr_act.max()
                        arr = arr_act / mx if mx > 0 else np.zeros(n)
            except Exception:
                pass
            return arr

        # Accumuler le score pondéré
        total_weight = sum(abs(c.get('weight', 0)) for c in config)
        if total_weight == 0:
            total_weight = 1.0

        score = np.zeros(n, dtype=float)
        for entry in config:
            mkey = entry.get('metric', '')
            w    = abs(entry.get('weight', 0)) / total_weight
            if w > 0:
                score += w * _get_layer(mkey)

        # Normalisation finale [0, 1]
        mx = score.max()
        if mx > 0:
            score /= mx

        mesh["Combined_Custom"] = score
        # Mémoriser la config pour l'affichage
        mesh.field_data['combined_config'] = np.array(
            [f"{c['metric']}*{c['weight']:.2f}" for c in config], dtype=object
        )

        # Stocker les temps d'activation sur le mesh (pour isochrones dans le renderer)
        for pacing, array_name in (('RV', 'Activation_RV_Time'), ('LV', 'Activation_LV_Time')):
            cache = f'activation_{pacing}'
            if cache in self.computed_metrics:
                src  = self.computed_metrics[cache]
                act  = src["Activation_Time"].astype(float)
                if len(act) == n:
                    mesh[array_name] = act
                else:
                    tree_src = cKDTree(src.points)
                    _, idx   = tree_src.query(mesh.points, k=1)
                    mesh[array_name] = act[idx]

        self.computed_metrics['combined_score'] = mesh
        return mesh

    def estimate_scar_from_thickness(self) -> None:
        """Estime la cicatrice depuis l'épaisseur (quand IRM manquant).
        Dense scar : T < 1 mm, border zone : 1–5 mm, total scar : T < 5 mm.
        L'estimation est grossière et sous-estime par rapport au LGE-MRI."""
        if 'wall_thickness' not in self.computed_metrics:
            return
        
        wt_mesh = self.computed_metrics['wall_thickness']
        T = wt_mesh["Wall_Thickness"].astype(float)
        
        # Créer masques basés sur seuils
        dense_mask = T < 1.0
        border_mask = (T >= 1.0) & (T <= 5.0)
        total_scar_mask = T < 5.0
        
        # Extraire régions comme pseudo-meshes
        # Dense scar
        if np.any(dense_mask):
            dense_mesh = wt_mesh.copy()
            dense_mesh = dense_mesh.threshold(value=1.0, scalars='Wall_Thickness', invert=True)
            self.meshes['dense_scar_estimated'] = dense_mesh
        
        # Border zone (via channels déjà calculé)
        if 'channels' in self.computed_metrics:
            border_mesh = self.computed_metrics['channels'].copy()
            self.meshes['border_zone_estimated'] = border_mesh
        
        # Total scar region
        if np.any(total_scar_mask):
            scar_mesh = wt_mesh.copy()
            scar_mesh = scar_mesh.threshold(value=5.0, scalars='Wall_Thickness', invert=True)
            self.meshes['scar_estimated'] = scar_mesh
        
        # Marquer que ce sont des estimations
        self._has_estimated_scar = True
    
    def _prepare_substrate(self, scar_decay: float = 5.0) -> pv.PolyData:
        """Prépare le substrat (CV continu) pour la simulation d'activation.

        La vitesse de conduction suit la formule sigmoïde de Cedilnik
        (EP-Europace 2018, éq. 3), purement basée sur l'épaisseur pariétale.
        Une décroissance exponentielle est appliquée aux abords de la cicatrice dense.

            v(X) = v_cedilnik(T) × exp(-d_to_dense_scar / scar_decay)

        où v_cedilnik(T) = v_max / (1 + exp(r * (p - t(X))))

        v_max = 0.6 m/s,  p = 3 mm (inflexion),  r = 2 (pente)
        
        Parameters
        ----------
        scar_decay : float
            Constante de décroissance (mm) pour exp(-d / decay). Défaut 5 mm.
            Values: 1–20 mm. Commandes des sliders /2 pour stabilité UI.
        """
        # Cache basé sur scar_decay
        cache_key = f'_substrate_decay_{scar_decay:.1f}'
        if cache_key in self.computed_metrics:
            return self.computed_metrics[cache_key]
        
        if 'wall_thickness' not in self.computed_metrics:
            self.compute_wall_thickness()
        
        mesh = self.computed_metrics['wall_thickness'].copy()
        pts  = mesh.points
        T    = mesh["Wall_Thickness"].astype(float)  # distance EPI (mm) = épaisseur pariétale
        n_pts = len(pts)
        
        # --- Formule sigmoïde continue de Cedilnik (EP-Europace 2018, éq. 3) ---
        #
        #   v(X) = v_max / (1 + exp(r * (p - t(X))))
        #
        #   v_max = 0.6 m/s   — vitesse du myocarde sain
        #   p     = 3.0 mm    — point d’inflexion (centre de la zone grise 2-5 mm)
        #   r     = 2         — pente : vitesse quasi nulle sous 2 mm (dense scar),
        #                                ~0.3 m/s à 3 mm (scar), ~0.56 m/s à 5 mm
        #
        # Tout est continu et dépend uniquement de l'épaisseur — pas de seuil discret.
        v_max = self.CV_HEALTHY   # 0.6 m/s
        p     = 3.0               # mm
        r     = 2.0               # dimensionless
        cv    = v_max / (1.0 + np.exp(r * (p - T)))
        
        # --- Décroissance exponentielle aux abords de dense_scar ---
        # Identifier les points dense_scar (T < 2 mm) et appliquer
        # un amortissement LOCALISÉ exp(-distance_to_dense_scar / scar_decay)
        # SEULEMENT pour les points à proximité (< 5 mm du dense scar)
        dense_mask = T < 2.0
        if np.any(dense_mask):
            dense_pts = pts[dense_mask]
            tree_dense = cKDTree(dense_pts)
            # Distance de chaque point au dense_scar le plus proche
            dists, _ = tree_dense.query(pts, k=1)
            # Modulation CV : appliquer SEULEMENT dans rayon 5 mm
            influence_radius = 5.0  # mm
            in_influence = dists <= influence_radius
            decay_mod = np.exp(-dists[in_influence] / (scar_decay + 1e-6))
            cv[in_influence] = cv[in_influence] * decay_mod
        
        # Stocker sur le mesh
        mesh["Conduction_Velocity_mm_ms"] = cv
        
        # --- Adjacence hybride pour Dijkstra éikonal ---
        #
        # Deux sources de voisins sont combinées :
        #   1. Arêtes du maillage : cohérence topologique
        #   2. Boule de 3 mm      : pont les bords ouverts laissés par le clipping
        #                           de l'anneau mitral/aortique, sans permettre de
        #                           franchir une zone cicatricielle (>3 mm d'étendue)
        #
        # Coût : dt = d / v_harm  (moyenne harmonique, pénalise fortement le slow tissue)
        
        # 1) Arêtes topologiques
        edge_pairs = set()
        raw_edges = mesh.extract_all_edges()
        for k in range(raw_edges.n_cells):
            cell = raw_edges.get_cell(k)
            if cell.n_points == 2:
                i, j = int(cell.point_ids[0]), int(cell.point_ids[1])
                edge_pairs.add((min(i, j), max(i, j)))
        
        # 2) Paires dans un rayon de 3 mm (comble les bords ouverts)
        tree_adj = cKDTree(pts)
        close_pairs = tree_adj.query_pairs(r=3.0, output_type='set')
        all_pairs = edge_pairs | close_pairs
        
        # Construire la liste d'adjacence
        adj = [[] for _ in range(len(pts))]
        for i, j in all_pairs:
            d = float(np.linalg.norm(pts[i] - pts[j]))
            if d < 1e-9:
                continue
            v_harm = 2.0 * cv[i] * cv[j] / (cv[i] + cv[j] + 1e-12)
            if v_harm < 1e-6:
                continue                    # arête bloquée (dense scar absolu)
            dt = d / v_harm
            adj[i].append((j, dt))
            adj[j].append((i, dt))
        
        mesh._adj = adj
        self.computed_metrics['_substrate'] = mesh
        return mesh
    
    def simulate_activation(self, pacing_site: str = "SR",
                            scar_decay: float = 5.0,
                            custom_point: np.ndarray = None) -> Optional[pv.PolyData]:
        """Simulation Eikonal (Dijkstra) depuis le site de pacing (SR/RV/LV/CUSTOM).
        Retourne le mesh avec Activation_Time et Conduction_Velocity.

        Parameters
        ----------
        scar_decay : float
            Constante de décroissance exponentielle du scar density (mm, défaut 5).
            Si différente de la dernière valeur, invalide le cache substrat.
        custom_point : np.ndarray, optional
            Coordonnées 3D du point de pacing choisi par l'utilisateur (clic droit).
            Si fourni, pacing_site est ignoré et ce point est utilisé.
        """
        substrate = self._prepare_substrate()
        mesh = substrate.copy()
        pts = mesh.points
        cv = mesh["Conduction_Velocity_mm_ms"]
        n_pts = len(pts)
        
        # Point de départ selon pacing_site (identique simu_conduct.py)
        if custom_point is not None:
            # User-selected pacing point: find nearest mesh vertex
            tree_custom = cKDTree(pts)
            _, start_idx = tree_custom.query(custom_point, k=1)
            start_idx = int(start_idx)
            pacing_site = "CUSTOM"
        elif pacing_site == "SR" or pacing_site == "SINUS":
            y_coords = pts[:, 1]
            x_coords = pts[:, 0]
            score = y_coords / np.std(y_coords) - x_coords / np.std(x_coords)
            start_idx = np.argmax(score)
        elif pacing_site == "RV":
            x_coords = pts[:, 0]
            z_coords = pts[:, 2]
            score = x_coords / np.std(x_coords) - np.abs(z_coords) / np.std(z_coords)
            start_idx = np.argmax(score)
        elif pacing_site == "LV":
            x_coords = pts[:, 0]
            y_coords = pts[:, 1]
            score = -x_coords / np.std(x_coords) - np.abs(y_coords - np.mean(y_coords)) / np.std(y_coords)
            start_idx = np.argmax(score)
        else:
            # Par défaut : point dans la zone la plus épaisse (tissu sain ≥ 5 mm)
            T_sub = substrate["Wall_Thickness"].astype(float)
            healthy = np.where(T_sub >= 5.0)[0]
            if len(healthy) > 0:
                start_idx = int(np.random.choice(healthy))
            else:
                start_idx = 0
        
        # Propagation Eikónal — Dijkstra sur arêtes de maillage
        # (pas de ball-tree : l'onde ne peut pas court-circuiter la cicatrice)
        adj = substrate._adj
        activation = np.full(n_pts, np.inf)
        activation[start_idx] = 0.0
        
        visited  = np.zeros(n_pts, dtype=bool)
        to_visit = [(0.0, start_idx)]
        heapq.heapify(to_visit)
        
        while to_visit:
            current_time, current = heapq.heappop(to_visit)
            
            if visited[current]:
                continue
            visited[current] = True
            
            for neighbor, dt in adj[current]:
                if visited[neighbor] or np.isinf(dt):
                    continue
                new_time = current_time + dt
                if new_time < activation[neighbor]:
                    activation[neighbor] = new_time
                    heapq.heappush(to_visit, (new_time, neighbor))
        
        # Gérer les points inaccessibles (identique simu_conduct.py)
        unreached = np.isinf(activation)
        if np.any(unreached) and np.any(~unreached):
            activation[unreached] = np.max(activation[~unreached]) + 50.0
        
        mesh["Activation_Time"] = activation
        mesh["Conduction_Velocity"] = cv

        # Stocker les coordonnées du point de pacing pour l'affichage
        self._pacing_points = getattr(self, '_pacing_points', {})
        self._pacing_points[pacing_site] = pts[start_idx].tolist()

        self.computed_metrics[f'activation_{pacing_site}'] = mesh

        return mesh

    def compute_local_entropy(self) -> Optional[pv.PolyData]:
        """Entropie locale de Shannon sur la distribution d'épaisseur.
        Utilise k-NN (k=30) + 8 bins pour un résultat lissé et cohérent."""
        if 'wall_thickness' not in self.computed_metrics:
            self.compute_wall_thickness()
        
        mesh = self.computed_metrics['wall_thickness'].copy()
        pts = mesh.points
        T = mesh["Wall_Thickness"].astype(float)
        n = len(T)
        
        # k-NN fixe pour éviter les artfacts de densité de maillage
        K = min(30, n - 1)
        tree = cKDTree(pts)
        _, knn_idx = tree.query(pts, k=K + 1)  # +1 car inclut le point lui-même
        
        N_BINS = 5
        t_min, t_max = np.min(T), np.max(T)
        if t_max - t_min < 1e-6:
            mesh["Local_Entropy"] = np.zeros(n)
            self.computed_metrics['local_entropy'] = mesh
            return mesh
        
        max_entropy = np.log2(N_BINS)
        
        entropy = np.zeros(n)
        for i in range(n):
            local_T = T[knn_idx[i]]  # k+1 voisins incluant le point
            counts, _ = np.histogram(local_T, bins=N_BINS, range=(t_min, t_max))
            probs = counts / counts.sum()
            probs = probs[probs > 0]
            entropy[i] = -np.sum(probs * np.log2(probs)) / max_entropy
        
        mesh["Local_Entropy"] = np.clip(entropy, 0, 1)
        self.computed_metrics['local_entropy'] = mesh
        return mesh
    
    # =====================================================================
    # Helpers pour le dashboard
    # =====================================================================
    def compute_all_metrics(self):
        """Pré-calcule toutes les métriques disponibles pour le dashboard"""
        if 'lv_epi_dist' in self.meshes:
            if 'wall_thickness' not in self.computed_metrics:
                self.compute_wall_thickness()
            if 'ciaccio_ratio' not in self.computed_metrics:
                self.compute_ciaccio_ratio()
            if 'deceleration_zones' not in self.computed_metrics:
                self.compute_deceleration_zones()
            if 'channels' not in self.computed_metrics:
                self.compute_channels()
            if 'isthmus' not in self.computed_metrics:
                self.compute_isthmus_prediction()
            if 'local_entropy' not in self.computed_metrics:
                self.compute_local_entropy()
            if 'laplacian' not in self.computed_metrics:
                self.compute_laplacian()
            if 'scar_burden' not in self.computed_metrics:
                self.compute_scar_burden()
            if 'channelness' not in self.computed_metrics:
                self.compute_channelness()
            if 'cv_map' not in self.computed_metrics:
                self.compute_cv_map()
            if 'border_zone_cedilnik' not in self.computed_metrics:
                self.compute_border_zone_cedilnik()
            if 'combined_zones' not in self.computed_metrics:
                self.compute_combined_zones()
        if 'scar_transmurality' in self.meshes:
            if 'transmurality' not in self.computed_metrics:
                self.compute_transmurality()

    def get_scalar_data(self, metric_key: str, scalar_name: str) -> Optional[np.ndarray]:
        """Retourne le tableau scalaire brut pour histogrammes/charts"""
        if metric_key in self.computed_metrics:
            mesh = self.computed_metrics[metric_key]
            if scalar_name in mesh.array_names:
                return mesh[scalar_name].astype(float).copy()
        return None

    # =====================================================================
    # Statistiques complètes
    # =====================================================================
    def get_statistics(self) -> Dict:
        """Calcule et retourne les statistiques globales complètes"""

        stats = {}

        # ═══════════════════════════════════════
        # GÉOMÉTRIE DU VENTRICULE
        # ═══════════════════════════════════════
        if 'lv_epi_dist' in self.meshes:
            mesh = self.meshes['lv_epi_dist']
            stats['LV_Surface_cm2'] = mesh.area / 100          # mm² → cm²
            stats['LV_N_Points']    = mesh.n_points
            stats['LV_N_Cells']     = mesh.n_cells

            bounds = mesh.bounds
            stats['LV_BBox_X_mm'] = bounds[1] - bounds[0]
            stats['LV_BBox_Y_mm'] = bounds[3] - bounds[2]
            stats['LV_BBox_Z_mm'] = bounds[5] - bounds[4]
            dims = [stats['LV_BBox_X_mm'], stats['LV_BBox_Y_mm'],
                    stats['LV_BBox_Z_mm']]
            stats['LV_Long_Axis_mm']  = max(dims)
            stats['LV_Short_Axis_mm'] = min(dims)

            # Volume (mesh fermé ou convex hull)
            try:
                vol = abs(mesh.volume)
                if vol > 0:
                    stats['LV_Volume_mL'] = vol / 1000         # mm³ → mL
            except Exception:
                pass
            if 'LV_Volume_mL' not in stats:
                try:
                    hull = ConvexHull(mesh.points)
                    stats['LV_Volume_mL'] = hull.volume / 1000
                except Exception:
                    pass

            # Indice de sphéricité
            if 'LV_Volume_mL' in stats and stats['LV_Volume_mL'] > 0:
                SA_mm2 = stats['LV_Surface_cm2'] * 100
                V_mm3  = stats['LV_Volume_mL'] * 1000
                r_eq   = (3 * V_mm3 / (4 * np.pi)) ** (1 / 3)
                SA_sph = 4 * np.pi * r_eq ** 2
                stats['LV_Sphericity'] = min(SA_sph / SA_mm2, 1.0)

        # ═══════════════════════════════════════
        # ÉPAISSEUR DE PAROI
        # ═══════════════════════════════════════
        if 'wall_thickness' in self.computed_metrics:
            mesh = self.computed_metrics['wall_thickness']
            T = mesh["Wall_Thickness"].astype(float)

            stats['T_mean']     = float(np.mean(T))
            stats['T_median']   = float(np.median(T))
            stats['T_std']      = float(np.std(T))
            stats['T_min']      = float(np.min(T))
            stats['T_max']      = float(np.max(T))
            stats['T_p5']       = float(np.percentile(T, 5))
            stats['T_p25']      = float(np.percentile(T, 25))
            stats['T_p75']      = float(np.percentile(T, 75))
            stats['T_p95']      = float(np.percentile(T, 95))
            stats['T_IQR']      = stats['T_p75'] - stats['T_p25']
            stats['T_CV']       = (stats['T_std'] / stats['T_mean']
                                   if stats['T_mean'] > 0 else 0)
            stats['T_skewness'] = float(scipy_skew(T))
            stats['T_kurtosis'] = float(scipy_kurtosis(T))

            # Classification tissulaire (par point)
            stats['Scar_Pct']        = float(np.mean(T < 1.0) * 100)
            stats['Border_Zone_Pct'] = float(np.mean((T >= 1.0) & (T <= 5.0)) * 100)
            stats['Healthy_Pct']     = float(np.mean(T > 5.0) * 100)
            stats['Thinning_Index']  = float(np.mean(T < 5.5) * 100)

            # Classification pondérée par la surface (plus précise)
            try:
                mesh_wt = self.computed_metrics['wall_thickness']
                cell_mesh  = mesh_wt.point_data_to_cell_data()
                cell_sizes = cell_mesh.compute_cell_sizes()
                cells_T    = cell_sizes.cell_data['Wall_Thickness']
                cell_areas = cell_sizes.cell_data['Area']
                total_area = float(np.sum(cell_areas))

                scar_a    = float(np.sum(cell_areas[cells_T < 1.0]))
                border_a  = float(np.sum(cell_areas[(cells_T >= 1.0)
                                                     & (cells_T <= 5.0)]))
                healthy_a = float(np.sum(cell_areas[cells_T > 5.0]))

                stats['Scar_Area_cm2']    = scar_a / 100
                stats['Border_Area_cm2']  = border_a / 100
                stats['Healthy_Area_cm2'] = healthy_a / 100
                if total_area > 0:
                    stats['Scar_Area_Pct']    = scar_a / total_area * 100
                    stats['Border_Area_Pct']  = border_a / total_area * 100
                    stats['Healthy_Area_Pct'] = healthy_a / total_area * 100
            except Exception:
                pass

        # ═══════════════════════════════════════
        # CIACCIO / RISQUE ARYTHMIQUE
        # ═══════════════════════════════════════
        if 'ciaccio_ratio' in self.computed_metrics:
            mesh = self.computed_metrics['ciaccio_ratio']
            rho = mesh["Ciaccio_Ratio"].astype(float)
            rho_clean = rho[rho < 5.0]

            if len(rho_clean) > 0:
                stats['Rho_Mean']   = float(np.mean(rho_clean))
                stats['Rho_Median'] = float(np.median(rho_clean))
                stats['Rho_Max']    = float(np.max(rho_clean))
                stats['Rho_Std']    = float(np.std(rho_clean))
                stats['Rho_p95']    = float(np.percentile(rho_clean, 95))

            stats['DZ_Area_Pct']  = float(np.mean(rho > stats.get('Rho_p95', 0.33)) * 100)

        # ═══════════════════════════════════════
        # ISTHMES
        # ═══════════════════════════════════════
        if 'isthmus' in self.computed_metrics:
            mesh = self.computed_metrics['isthmus']
            isth = mesh["Isthmus_Probability"]
            stats['Isthmus_N_Points'] = int(np.sum(isth > 0.5))
            stats['Isthmus_Pct']      = float(np.mean(isth > 0.5) * 100)
            if "Isthmus_Boundaries" in mesh.array_names:
                stats['Isthmus_Boundary_Points'] = int(
                    np.sum(mesh["Isthmus_Boundaries"] > 0.5))

        # ═══════════════════════════════════════
        # CANAUX DE CONDUCTION
        # ═══════════════════════════════════════
        if 'channels' in self.computed_metrics:
            mesh = self.computed_metrics['channels']
            ch = mesh["Channel_Zone"]
            stats['Channel_Pct'] = float(np.mean(ch > 0.5) * 100)
            cs = mesh["Channel_Score"]
            cs_active = cs[cs > 0]
            if len(cs_active) > 0:
                stats['Channel_Score_Mean'] = float(np.mean(cs_active))
                stats['Channel_Score_Max']  = float(np.max(cs_active))

        # ═══════════════════════════════════════
        # SCAR BURDEN (Utah) — depuis épaisseur
        # ═══════════════════════════════════════
        if hasattr(self, '_scar_burden') and self._scar_burden:
            sb = self._scar_burden
            stats['Scar_Burden_Pct']    = sb.get('utah_pct',    sb.get('scar_burden_pct', 0.0))
            stats['Dense_Burden_Pct']   = sb.get('dense_pct',   sb.get('dense_burden_pct', 0.0))
            stats['Utah_Grade']         = sb.get('utah_label',  '')
            stats['Healthy_Area_CT_pct']= sb.get('healthy_pct', 0.0)
            stats['Border_Area_CT_pct'] = sb.get('border_pct',  0.0)
            stats['Scar_Area_CT_pct']   = sb.get('scar_pct',    0.0)
            stats['Dense_Area_CT_pct']  = sb.get('dense_pct',   0.0)
            stats['Total_Area_CT_cm2']  = sb.get('total_area_cm2', 0.0)
        
        # ═══════════════════════════════════════
        # CHANNELNESS
        # ═══════════════════════════════════════
        if 'channelness' in self.computed_metrics:
            mesh_ch = self.computed_metrics['channelness']
            ch_vals = mesh_ch["Channelness"].astype(float)
            ch_active = ch_vals[ch_vals > 0]
            
            if len(ch_active) > 0:
                stats['Channelness_Mean'] = float(np.mean(ch_active))
                stats['Channelness_Max'] = float(np.max(ch_active))
                stats['Channelness_Median'] = float(np.median(ch_active))
                stats['Channelness_Std'] = float(np.std(ch_active))
                stats['Channelness_High_Pct'] = float(np.mean(ch_vals > 0.5) * 100)
                stats['Channelness_P95'] = float(np.percentile(ch_active, 95))
            
            if "Channel_Width_mm" in mesh_ch.array_names:
                widths = mesh_ch["Channel_Width_mm"]
                widths_active = widths[widths > 0]
                if len(widths_active) > 0:
                    stats['Channel_Width_Mean_mm'] = float(np.mean(widths_active))
                    stats['Channel_Width_Min_mm'] = float(np.min(widths_active))
        
        # Anatomical channelness statistics
        if 'anatomical_channelness' in self.computed_metrics:
            mesh_ac = self.computed_metrics['anatomical_channelness']
            if 'Anatomical_Channelness' in mesh_ac.array_names:
                ac_vals = mesh_ac['Anatomical_Channelness'].astype(float)
                ac_active = ac_vals[ac_vals > 0.05]
                if len(ac_active) > 0:
                    stats['AnatomCh_Mean']    = float(np.mean(ac_active))
                    stats['AnatomCh_Max']     = float(np.max(ac_active))
                    stats['AnatomCh_High_Pct']= float(np.mean(ac_vals > 0.5) * 100)
            if 'narrow_pct' in mesh_ac.field_data:
                stats['AnatomCh_Narrow_Pct'] = float(mesh_ac.field_data['narrow_pct'][0])
            if 'max_lw_raw' in mesh_ac.field_data:
                stats['AnatomCh_MaxLW']      = float(mesh_ac.field_data['max_lw_raw'][0])
            if 'n_skeleton_pts' in mesh_ac.field_data:
                stats['AnatomCh_SkelPts']    = int(mesh_ac.field_data['n_skeleton_pts'][0])

        if hasattr(self, '_cedilnik_stats') and self._cedilnik_stats:
            cs = self._cedilnik_stats
            stats['Cedilnik_BZ_Pct'] = cs.get('bz_point_pct', 0.0)
            stats['Cedilnik_BZ_Area_cm2'] = cs.get('bz_area_cm2', 0.0)
            stats['Cedilnik_BZ_Area_Pct'] = cs.get('bz_area_pct', 0.0)
            stats['Cedilnik_Radius_mm'] = cs.get('dilation_radius_mm', 2.0)
        
        # ═══════════════════════════════════════
        # CV MAP (Cedilnik transfer function)
        # ═══════════════════════════════════════
        if hasattr(self, '_cv_stats') and self._cv_stats:
            cvs = self._cv_stats
            stats['CV_Mean_ms'] = cvs.get('cv_mean', 0.0)
            stats['CV_Min_ms'] = cvs.get('cv_min', 0.0)
            stats['CV_Max_ms'] = cvs.get('cv_max', 0.0)
            stats['CV_Blocked_Pct'] = cvs.get('pct_blocked', 0.0)
            stats['CV_Slow_Pct'] = cvs.get('pct_slow', 0.0)
            stats['CV_Normal_Pct'] = cvs.get('pct_normal', 0.0)
        
        # ═══════════════════════════════════════
        # ═══════════════════════════════════════
        # COMBINED ZONES
        # ═══════════════════════════════════════
        if 'combined_zones' in self.computed_metrics:
            mesh_cz = self.computed_metrics['combined_zones']
            cz = mesh_cz["Combined_Score"].astype(float)
            stats['Combined_Score_Mean'] = float(np.mean(cz))
            stats['Combined_High_Risk_Pct'] = float(np.mean(cz >= 0.6) * 100)
            stats['Combined_Critical_Pct'] = float(np.mean(cz >= 0.8) * 100)

        # ═══════════════════════════════════════
        # CICATRICE — IRM ou ESTIMATION
        # ═══════════════════════════════════════
        # Prioriser IRM si disponible, sinon utiliser estimation depuis T
        scar_source = None  # 'mri' ou 'estimated'
        
        if 'dense_scar' in self.meshes:
            scar_mesh = self.meshes['dense_scar']
            scar_source = 'mri'
        elif 'dense_scar_estimated' in self.meshes:
            scar_mesh = self.meshes['dense_scar_estimated']
            scar_source = 'estimated'
        
        if scar_source:
            stats['Dense_Scar_Surface_cm2'] = scar_mesh.area / 100
            stats['Dense_Scar_N_Points']    = scar_mesh.n_points
            stats['Scar_Data_Source']       = scar_source  # Pour traçabilité

            # Scar burden = surface cicatrice / surface VG totale
            if 'LV_Surface_cm2' in stats and stats['LV_Surface_cm2'] > 0:
                stats['Dense_Scar_Burden_Pct'] = (
                    stats['Dense_Scar_Surface_cm2']
                    / stats['LV_Surface_cm2'] * 100)

            # Volume tissulaire = intégrale(épaisseur × dA) sur la surface
            if 'wall_thickness' in self.computed_metrics:
                wt_mesh = self.computed_metrics['wall_thickness']
                T_all   = wt_mesh["Wall_Thickness"].astype(float)
                tree_wt = cKDTree(wt_mesh.points)

                try:
                    cc = scar_mesh.cell_centers()
                    _, idx = tree_wt.query(cc.points, k=1)
                    T_at_scar_cells = T_all[idx]

                    scar_sizes = scar_mesh.compute_cell_sizes()
                    scar_areas = scar_sizes.cell_data['Area']

                    tissue_vol = float(np.sum(scar_areas * T_at_scar_cells))
                    stats['Dense_Scar_Tissue_Volume_cm3'] = tissue_vol / 1000

                    # Épaisseur moyenne dans la cicatrice
                    _, idx_pts = tree_wt.query(scar_mesh.points, k=1)
                    stats['Scar_Mean_Thickness_mm'] = float(
                        np.mean(T_all[idx_pts]))
                except Exception:
                    pass

            # Compacité de la cicatrice
            try:
                edges = scar_mesh.extract_feature_edges(
                    boundary_edges=True, manifold_edges=False,
                    feature_edges=False, non_manifold_edges=False)
                if edges.n_points > 0 and hasattr(edges, 'length'):
                    perimeter = edges.length
                    if perimeter > 0:
                        stats['Scar_Compactness'] = (
                            4 * np.pi * scar_mesh.area / perimeter ** 2)
            except Exception:
                pass

        # Scar (LE) total - IRM ou estimation
        scar_le_source = None
        if 'scar' in self.meshes:
            scar_le = self.meshes['scar']
            scar_le_source = 'mri'
        elif 'scar_estimated' in self.meshes:
            scar_le = self.meshes['scar_estimated']
            scar_le_source = 'estimated'
        
        if scar_le_source:
            stats['Scar_LE_Surface_cm2'] = scar_le.area / 100
            if 'LV_Surface_cm2' in stats and stats['LV_Surface_cm2'] > 0:
                stats['Scar_LE_Burden_Pct'] = (
                    stats['Scar_LE_Surface_cm2']
                    / stats['LV_Surface_cm2'] * 100)
            if scar_le_source == 'estimated':
                stats['Scar_LE_Data_Source'] = 'estimated'

        # Distribution endo/intra/epi (uniquement si IRM réelle)
        for key, label in [('scar_endo', 'Endo'),
                           ('scar_intra', 'Intra'),
                           ('scar_epi', 'Epi')]:
            if key in self.meshes:
                stats[f'Scar_{label}_Surface_cm2'] = (
                    self.meshes[key].area / 100)

        # ═══════════════════════════════════════
        # TRANSMURALITÉ
        # ═══════════════════════════════════════
        if 'transmurality' in self.computed_metrics:
            mesh  = self.computed_metrics['transmurality']
            trans = mesh["Transmurality"].astype(float)

            stats['Trans_Mean']   = float(np.mean(trans))
            stats['Trans_Median'] = float(np.median(trans))
            stats['Trans_Std']    = float(np.std(trans))
            stats['Trans_Min']    = float(np.min(trans))
            stats['Trans_Max']    = float(np.max(trans))
            stats['Trans_p25']    = float(np.percentile(trans, 25))
            stats['Trans_p75']    = float(np.percentile(trans, 75))

            trans_nz = trans[trans > 0]
            if len(trans_nz) > 0:
                stats['Trans_NZ_Mean']   = float(np.mean(trans_nz))
                stats['Trans_NZ_Median'] = float(np.median(trans_nz))

            total = len(trans)
            stats['Trans_None_Pct']      = float(
                np.sum(trans <= 0) / total * 100)
            stats['Trans_Subendo_Pct']   = float(
                np.sum((trans > 0) & (trans <= 25)) / total * 100)
            stats['Trans_Midmural_Pct']  = float(
                np.sum((trans > 25) & (trans <= 50)) / total * 100)
            stats['Trans_Subepi_Pct']    = float(
                np.sum((trans > 50) & (trans <= 75)) / total * 100)
            stats['Trans_Transmural_Pct'] = float(
                np.sum(trans > 75) / total * 100)

        # ═══════════════════════════════════════
        # ENTROPIE LOCALE
        # ═══════════════════════════════════════
        if 'local_entropy' in self.computed_metrics:
            mesh = self.computed_metrics['local_entropy']
            ent  = mesh["Local_Entropy"].astype(float)
            stats['Entropy_Mean']        = float(np.mean(ent))
            stats['Entropy_Median']      = float(np.median(ent))
            stats['Entropy_Max']         = float(np.max(ent))
            stats['Entropy_Std']         = float(np.std(ent))
            # P75 et P95 plutôt qu'un seuil arbitraire
            stats['Entropy_P75']         = float(np.percentile(ent, 75))
            stats['Entropy_P95']         = float(np.percentile(ent, 95))

        # ═══════════════════════════════════════
        # GRAISSE
        # ═══════════════════════════════════════
        if 'lv_fat' in self.meshes:
            fat_mesh = self.meshes['lv_fat']
            stats['Fat_Surface_cm2'] = fat_mesh.area / 100
            if 'LV_Surface_cm2' in stats and stats['LV_Surface_cm2'] > 0:
                stats['Fat_Burden_Pct'] = (
                    stats['Fat_Surface_cm2']
                    / stats['LV_Surface_cm2'] * 100)

        # Evidence-based risk score (continuous metrics, literature thresholds)
        risk = {}

        # 1. Scar Burden (threshold ≥10%, Ponnusamy 2023)
        scar_burden = stats.get('Scar_Burden_Pct',
                                stats.get('Dense_Scar_Burden_Pct',
                                stats.get('Scar_Area_Pct', 0.0)))
        risk['scar_burden_pct'] = scar_burden
        risk['scar_burden_norm'] = min(scar_burden / 20.0, 1.0)  # normalisé 0-1
        risk['scar_burden_ref'] = "Ponnusamy 2023 (PMID 37217065): ≥10%"

        # 2. Scar Volume (threshold ≥37.3 mL, John 2023 PAINES2D)
        scar_vol = stats.get('Dense_Scar_Tissue_Volume_cm3', 0.0)
        risk['scar_volume_cm3'] = scar_vol
        risk['scar_volume_norm'] = min(scar_vol / 37.3, 1.0)  # normalisé sur seuil
        risk['scar_volume_ref'] = "John 2023 (PMID 37354175): ≥37.3 mL"

        # 3. Wall Thinning (proxy: % surface <1mm, Marchlinski 2000)
        thin_pct = stats.get('Scar_Area_Pct', stats.get('Scar_Pct', 0.0))
        risk['thinning_pct'] = thin_pct
        risk['thinning_norm'] = min(thin_pct / 15.0, 1.0)
        risk['thinning_ref'] = "Marchlinski 2000 (PMID 10725289): proxy épaisseur"

        # 4. Deceleration Zones (Raiman & Tung 2018)
        dz_pct = stats.get('DZ_Area_Pct', 0.0)
        risk['dz_extent_pct'] = dz_pct
        risk['dz_extent_norm'] = min(dz_pct / 20.0, 1.0)
        risk['dz_ref'] = "Raiman & Tung 2018 (PMID 30033360): CV < 0.6 m/s"

        # 5. Channels (exploratory, no validated threshold)
        ch_pct = stats.get('Channel_Pct', 0.0)
        risk['channel_extent_pct'] = ch_pct
        risk['channel_extent_norm'] = min(ch_pct / 30.0, 1.0)
        risk['channel_ref'] = "Pas de seuil validé — variable exploratoire"

        # 6. Isthmus (exploratory)
        isth_pct = stats.get('Isthmus_Pct', 0.0)
        risk['isthmus_pct'] = isth_pct
        risk['isthmus_norm'] = min(isth_pct / 10.0, 1.0)
        risk['isthmus_ref'] = "Pas de seuil validé — variable exploratoire"

        # 7. Transmurality (exploratory, ECG-guided decision)
        trans_pct = stats.get('Trans_Transmural_Pct', 0.0)
        risk['transmural_pct'] = trans_pct
        risk['transmural_norm'] = min(trans_pct / 30.0, 1.0)
        risk['transmural_ref'] = "Pas de seuil validé — décision clinique ECG-guided"

        # 8. Entropy (exploratory)
        ent_mean = stats.get('Entropy_Mean', 0.0)
        risk['entropy_mean'] = ent_mean
        risk['entropy_norm'] = min(ent_mean / 0.8, 1.0)
        risk['entropy_ref'] = "Pas de seuil validé — variable exploratoire"

        # 9. Sphericity Index (threshold ≥0.70)
        sph = stats.get('LV_Sphericity', 0.0)
        risk['sphericity'] = sph
        risk['sphericity_norm'] = min(sph / 0.70, 1.0)  # normalisé sur seuil 0.70
        risk['sphericity_ref'] = "Seuil ≥0.70 (remodelage ventriculaire)"

        # 10. Channelness (Cedilnik, exploratory)
        ch_high = stats.get('Channelness_High_Pct', 0.0)
        risk['channelness_high_pct'] = ch_high
        risk['channelness_norm'] = min(ch_high / 15.0, 1.0)
        risk['channelness_ref'] = "Cedilnik, Europace 2018 — variable exploratoire"

        # 11. Border Zone Cedilnik (exploratory)
        bz_pct = stats.get('Cedilnik_BZ_Area_Pct', 0.0)
        risk['bz_cedilnik_pct'] = bz_pct
        risk['bz_cedilnik_norm'] = min(bz_pct / 30.0, 1.0)
        risk['bz_cedilnik_ref'] = "Cedilnik, JACC:EP 2023 — variable exploratoire"

        # Composite score (weighted sum, EXPLORATORY — not clinically validated)
        weights = {
            'scar_burden':  3,  # littérature : seuil ≥10%
            'scar_volume':  3,  # littérature : seuil ≥37.3 mL
            'thinning':     2,  # corrélation thickness–voltage établie
            'dz_extent':    3,  # littérature : DZ = sites ablation réussis
            'channel':      1,  # exploratoire
            'isthmus':      1,  # exploratoire
            'transmural':   2,  # corrélation avec approche (endo/épi)
            'entropy':      1,  # exploratoire
            'sphericity':   2,  # seuil 0.70, remodelage ventriculaire
            'channelness':  1,  # exploratoire (Cedilnik)
            'bz_cedilnik':  1,  # exploratoire (Cedilnik)
        }
        norms = [
            risk['scar_burden_norm'],
            risk['scar_volume_norm'],
            risk['thinning_norm'],
            risk['dz_extent_norm'],
            risk['channel_extent_norm'],
            risk['isthmus_norm'],
            risk['transmural_norm'],
            risk['entropy_norm'],
            risk['sphericity_norm'],
            risk['channelness_norm'],
            risk['bz_cedilnik_norm'],
        ]
        w_list = list(weights.values())
        total_weight = sum(w_list)
        composite = sum(n * w for n, w in zip(norms, w_list)) / total_weight

        risk['composite_score'] = composite  # 0-1
        risk['composite_level'] = (
            'Low' if composite < 0.25 else
            'Moderate' if composite < 0.50 else
            'High' if composite < 0.75 else
            'Very High'
        )
        risk['composite_note'] = (
            "EXPLORATOIRE — Non validé cliniquement. "
            "Pondérations basées sur le niveau de preuve de chaque composante."
        )

        stats['_risk'] = risk

        return stats