"""
Gestionnaire de visualisation - Gère les rendus PyVista
Fidèle aux rendus de LV_topology_v.0.2.py, simu_conduct.py, advanced_scar_analysis.py
"""

import pyvista as pv
import numpy as np
from typing import Optional, Dict
from scipy.spatial import cKDTree
from ..metrics.metrics_catalog import MetricsCatalog, MetricDefinition


class VisualizationManager:
    """Gère le rendu des métriques avec PyVista, support multi-vues simultanées"""
    
    def __init__(self, plotter):
        self.plotter = plotter
        self.catalog = MetricsCatalog()
        self.current_meshes = {}       # Cache des meshes par vue
        self.view_actors = {}          # Acteurs par vue pour nettoyage ciblé
        self._base_mesh = None         # Ghost mesh EPI pour overlays
        self._cursor_enabled = False   # Curseur synchronisé
        self._cursor_actors = {}       # Acteurs sphère curseur par vue key
        self._mesh_trees = {}          # KDTree par vue key pour snap rapide
        self._cursor_obs_id = None     # ID de l'observer VTK
        
    def _view_key(self, subplot_index: tuple) -> str:
        return f"{subplot_index[0]}_{subplot_index[1]}"

    # ─────────────────────────────────────────
    # CURSEUR SYNCHRONISÉ
    # ─────────────────────────────────────────
    def enable_linked_cursor(self):
        """Active le curseur synchronisé sur toutes les vues.
        Quand la souris survole n'importe quelle vue, une sphère lumineuse
        apparaît à la position correspondante sur chaque autre vue."""
        if self._cursor_enabled:
            return
        try:
            iren = self.plotter.iren
            vtk_iren = iren.interactor if hasattr(iren, 'interactor') else iren
            self._cursor_obs_id = vtk_iren.AddObserver(
                'MouseMoveEvent', self._on_mouse_move
            )
            self._cursor_enabled = True
        except Exception:
            pass

    def disable_linked_cursor(self):
        """Désactive le curseur synchronisé."""
        if not self._cursor_enabled:
            return
        try:
            iren = self.plotter.iren
            vtk_iren = iren.interactor if hasattr(iren, 'interactor') else iren
            if self._cursor_obs_id is not None:
                vtk_iren.RemoveObserver(self._cursor_obs_id)
            self._cursor_obs_id = None
        except Exception:
            pass
        self._remove_all_cursor_actors()
        self._cursor_enabled = False

    def _rebuild_mesh_trees(self):
        """Reconstruit les KDTrees pour tous les meshes actuellement affichés."""
        self._mesh_trees = {}
        for key, mesh in self.current_meshes.items():
            if mesh is not None and mesh.n_points > 0:
                try:
                    self._mesh_trees[key] = cKDTree(mesh.points)
                except Exception:
                    pass

    def _on_mouse_move(self, vtk_iren, event_name):
        """Callback VTK : pick la position 3D et place les sphères curseur."""
        if not self.current_meshes:
            return
        try:
            # Coordonnées écran
            x, y = vtk_iren.GetEventPosition()
            rw = vtk_iren.GetRenderWindow()

            # Picker sur toute la fenêtre
            import vtk
            picker = vtk.vtkWorldPointPicker()
            picker.Pick(x, y, 0, rw.GetRenderers().GetFirstRenderer())
            pos = np.array(picker.GetPickPosition())

            # Si le picker retourne l'origine (rien touché), on cache
            if np.allclose(pos, 0.0, atol=1e-6):
                self._remove_all_cursor_actors()
                return

            # Pour chaque vue, snapper la position au point le plus proche du mesh
            rows, cols = self.plotter.shape
            for r in range(rows):
                for c in range(cols):
                    key = f'{r}_{c}'
                    mesh = self.current_meshes.get(key)
                    if mesh is None:
                        continue
                    # S'assurer que le KDTree existe
                    if key not in self._mesh_trees:
                        try:
                            self._mesh_trees[key] = cKDTree(mesh.points)
                        except Exception:
                            continue
                    # Trouver le point le plus proche
                    _, idx = self._mesh_trees[key].query(pos, k=1)
                    snap_pos = mesh.points[idx]
                    self._draw_cursor(key, snap_pos, (r, c))
        except Exception:
            pass

    def _draw_cursor(self, key: str, pos: np.ndarray, subplot_index: tuple):
        """Place ou déplace la sphère curseur dans le subplot donné."""
        actor_name = f'_cursor_{key}'
        try:
            sphere = pv.Sphere(radius=1.5, center=pos)
            self.plotter.subplot(*subplot_index)
            self.plotter.add_mesh(
                sphere,
                name=actor_name,
                color='cyan',
                ambient=0.8,
                diffuse=0.3,
                specular=1.0,
                specular_power=60,
                render_points_as_spheres=False,
                show_scalar_bar=False,
                reset_camera=False,
            )
            self._cursor_actors[key] = actor_name
        except Exception:
            pass

    def _remove_all_cursor_actors(self):
        """Supprime toutes les sphères curseur."""
        for actor_name in self._cursor_actors.values():
            try:
                self.plotter.remove_actor(actor_name)
            except Exception:
                pass
        self._cursor_actors = {}
    
    def _clear_subplot(self, subplot_index: tuple):
        """Nettoie UNIQUEMENT les acteurs d'un subplot donné, sans toucher aux autres"""
        i, j = subplot_index
        key = self._view_key(subplot_index)
        
        self.plotter.subplot(i, j)
        
        # Supprimer les acteurs nommés de cette vue
        if key in self.view_actors:
            for actor_name in self.view_actors[key]:
                try:
                    self.plotter.remove_actor(actor_name)
                except:
                    pass
            self.view_actors[key] = []
        
        # Aussi supprimer les acteurs génériques de ce renderer
        renderer = self.plotter.renderers[i * self.plotter.shape[1] + j]
        renderer.RemoveAllViewProps()
        
    def _add_tracked_mesh(self, subplot_index: tuple, mesh, name_suffix: str, **kwargs):
        """Ajoute un mesh avec un nom unique par subplot pour tracking"""
        key = self._view_key(subplot_index)
        actor_name = f"{key}_{name_suffix}"
        
        if key not in self.view_actors:
            self.view_actors[key] = []
        self.view_actors[key].append(actor_name)
        
        self.plotter.add_mesh(mesh, name=actor_name, **kwargs)
    
    def set_clim_on_actors(self, subplot_index: tuple, clim: tuple):
        """Met à jour le scalar_range (clim) de tous les acteurs du subplot sans re-rendre.
        Utile pour les changements de slider qui doivent être fluides."""
        i, j = subplot_index
        try:
            self.plotter.subplot(i, j)
            # Récupérer le renderer du subplot courant
            renderer = self.plotter.renderer
            
            # Mettre à jour tous les acteurs du renderer
            for actor in renderer.GetActors():
                if actor is not None and hasattr(actor, 'mapper') and actor.mapper is not None:
                    if hasattr(actor.mapper, 'SetScalarRange'):
                        actor.mapper.SetScalarRange(*clim)
        except Exception:
            pass
        
    def render_metric(self, subplot_index: tuple, mesh: pv.PolyData, 
                     metric_def: MetricDefinition, dz_threshold: float = 0.33,
                     lap_threshold: float = 0.10, lap_min_area: float = 0.0):
        """
        Rend une métrique dans un subplot spécifique.
        Utilise les techniques de rendu des scripts de référence:
        - threshold() pour DZ, channels, isthmus
        - ghost mesh + overlay coloré
        - scalar_bar avec paramètres fidèles
        """
        i, j = subplot_index
        key = self._view_key(subplot_index)
        
        # 1. Nettoyer UNIQUEMENT ce subplot
        self._clear_subplot(subplot_index)
        self.plotter.subplot(i, j)
        
        # Sauvegarder le mesh et mettre à jour le KDTree pour le curseur
        self.current_meshes[key] = mesh
        if self._cursor_enabled and mesh is not None and mesh.n_points > 0:
            try:
                self._mesh_trees[key] = cKDTree(mesh.points)
            except Exception:
                pass
        
        # 2. Rendu selon le type de métrique (fidèle aux scripts de référence)
        metric_id = metric_def.id
        
        if metric_id == 'thickness' or metric_id == 'parietal_thickness':
            self._render_thickness(subplot_index, mesh, metric_def)
        elif metric_id == 'ciaccio':
            self._render_ciaccio(subplot_index, mesh, metric_def)
        elif metric_id == 'deceleration':
            self._render_deceleration_zones(subplot_index, mesh, metric_def, dz_threshold=dz_threshold)
        elif metric_id == 'channels':
            self._render_channels(subplot_index, mesh, metric_def)
        elif metric_id == 'laplacian':
            self._render_laplacian(subplot_index, mesh, metric_def, lap_threshold=lap_threshold, min_area_mm2=lap_min_area)
        elif metric_id == 'isthmus':
            self._render_isthmus(subplot_index, mesh, metric_def)
        elif metric_id == 'channelness':
            self._render_channelness(subplot_index, mesh, metric_def)
        elif metric_id == 'anatomical_channelness':
            self._render_anatomical_channelness(subplot_index, mesh, metric_def)
        elif metric_id == 'transmurality':
            self._render_transmurality(subplot_index, mesh, metric_def)
        elif metric_id in ('activation_sr', 'activation_rv', 'activation_lv', 'simulation'):
            self._render_activation(subplot_index, mesh, metric_def)
        elif metric_id == 'local_entropy':
            self._render_local_entropy(subplot_index, mesh, metric_def)
        elif metric_id == 'tri':
            self._render_tri(subplot_index, mesh, metric_def)
        elif metric_id == 'scar_distribution':
            self._render_scar_distribution(subplot_index, mesh, metric_def)
        elif metric_id == 'scar_burden':
            self._render_scar_burden(subplot_index, mesh, metric_def)
        elif metric_id == 'cv_map':
            self._render_cv_map(subplot_index, mesh, metric_def)
        elif metric_id == 'combined_score':
            self._render_combined_score(subplot_index, mesh, metric_def)
        else:
            # Fallback générique
            self._render_generic(subplot_index, mesh, metric_def)
        
        # 3. Titre (comme les scripts de référence, couleur fidèle)
        self.plotter.add_text(
            metric_def.name,
            position='upper_left',
            font_size=10,
            color=metric_def.title_color,
            name=f"{key}_title"
        )
    
    # =========================================================================
    # RENDUS TOPOGRAPHIE — identiques à LV_topology_v.0.2.py
    # =========================================================================
    
    def _render_thickness(self, si, mesh, mdef):
        """Rendu d'épaisseur de paroi avec 6 couleurs discrètes.
        Le clim_max est fourni par le slider (mdef.clim[1]), défaut 6 mm.
        Supporte 'thickness' (EPI_Distance) et 'parietal_thickness' (Wall_Thickness)."""
        self.plotter.subplot(*si)
        import numpy as np
        from matplotlib import cm
        from matplotlib.colors import ListedColormap

        # Lire clim_max depuis la définition (pilotée par le slider)
        clim_max = float(mdef.clim[1]) if (mdef.clim is not None and mdef.clim[1] > 0) else 6.0

        # 6 couleurs — colormap 'hot' échantillonnée de 0.03 à 0.98 (évite le blanc pur)
        n = 6
        hot = cm.get_cmap('hot', 256)
        colors = [hot(i / (n - 1) * 0.95 + 0.03) for i in range(n)]
        discrete_cmap = ListedColormap(colors)

        step = clim_max / n
        sbar = {
            'title': f'mm  (\u00d7{step:.2f}mm)',
            'color': 'black',
            'n_labels': n + 1,
            'vertical': True,
            'title_font_size': 10,
            'label_font_size': 8,
        }
        scalar_name = mdef.scalar_name  # EPI_Distance ou Wall_Thickness
        self._add_tracked_mesh(si, mesh, mdef.id,
            scalars=scalar_name,
            cmap=discrete_cmap,
            clim=[0, clim_max],
            scalar_bar_args=sbar,
            show_scalar_bar=True)
    
    def _render_ciaccio(self, si, mesh, mdef):
        """Vue 2 de LV_topology — IDENTIQUE :
        p.add_mesh(mesh_gradient, scalars="Ciaccio_Ratio_Display", cmap="jet_r", clim=[0, 0.5])
        """
        self.plotter.subplot(*si)
        sbar = {
            'title': 'Wavefront Curvature',
            'color': 'black',
            'n_labels': 5,
            'vertical': True,
            'title_font_size': 10,
            'label_font_size': 8
        }
        self._add_tracked_mesh(si, mesh, "ciaccio",
            scalars="Ciaccio_Ratio_Display",
            cmap="jet_r",
            clim=[0, 0.5],
            scalar_bar_args=sbar,
            show_scalar_bar=True)
    
    def _render_local_entropy(self, si, mesh, mdef):
        """Entropie locale — identification des zones chaotiques/irrégulières
        Échelle discrétisée en 5 niveaux."""
        self.plotter.subplot(*si)
        from matplotlib.colors import ListedColormap
        # Discretized entropy: 5 levels
        colors = [
            '#1a0533',  # Très faible (sombre)
            '#721f81',  # Faible
            '#f1605d',  # Modéré
            '#feb078',  # Élevé
            '#fcfdbf',  # Très élevé (clair)
        ]
        discrete_cmap = ListedColormap(colors)
        sbar = {'title': 'Entropy', 'color': 'black', 'n_labels': 5,
                'vertical': True, 'title_font_size': 10, 'label_font_size': 8}
        self._add_tracked_mesh(si, mesh, "local_entropy",
            scalars="Local_Entropy",
            cmap=discrete_cmap,
            clim=[0, 1],
            scalar_bar_args=sbar,
            show_scalar_bar=True)
    
    def _render_tri(self, si, mesh, mdef):
        """TRI (Terrain Ruggedness Index) — rugosité de surface basée sur variations d'épaisseur
        5 classes (suppression de la classe extrême blanche)."""
        self.plotter.subplot(*si)
        
        tri = mesh["TRI"].astype(float)
        tri_p95 = float(np.percentile(tri, 95)) if np.any(tri > 0) else 1.0
        
        # Colormap discrète qualitative — 5 classes
        from matplotlib.colors import ListedColormap
        colors = [
            '#2d7f5e',  # Vert foncé (très lisse, TRI faible)
            '#5dade2',  # Bleu clair (lisse)
            '#f39c12',  # Orange (modéré)
            '#e74c3c',  # Rouge (rugueux)
            '#8e44ad',  # Violet (très rugueux)
        ]
        discrete_cmap = ListedColormap(colors)
        
        sbar = {'title': 'TRI (mm)', 'color': 'black', 'n_labels': 5,
                'vertical': True, 'title_font_size': 10, 'label_font_size': 8}
        
        self._add_tracked_mesh(si, mesh, "tri",
            scalars="TRI",
            cmap=discrete_cmap,
            clim=[0, tri_p95],
            scalar_bar_args=sbar,
            show_scalar_bar=True)
    
    def _render_deceleration_zones(self, si, mesh, mdef, dz_threshold=0.33):
        """Vue 3 de LV_topology :
        Ghost mesh gris + zones de décélération en rouge, seuil paramétrable.
        """
        self.plotter.subplot(*si)
        
        # Ghost mesh (base grise transparente)
        self._add_tracked_mesh(si, mesh, "dz_ghost",
            color="lightgrey", opacity=0.1, show_scalar_bar=False)
        
        # Zones de décélération : threshold sur Ciaccio_Ratio
        try:
            dz = mesh.threshold(dz_threshold, scalars="Ciaccio_Ratio")
            if dz.n_points > 0:
                self._add_tracked_mesh(si, dz, "dz_zones",
                    color="red", show_scalar_bar=False)
        except:
            pass
    
    def _render_channels(self, si, mesh, mdef):
        """Canaux de conduction :
        UN SEUL MESH avec colormap catégorielle → élimine le z-fighting.
          - Canaux 0-5mm (colormap symétrique bleu→rouge@2.5mm→bleu)
          - Zone bordante (orange)
          - Reste (gris clair)
        """
        self.plotter.subplot(*si)
        
        # Si l'utilisateur a demandé une mise en évidence d'épaisseur cible,
        # afficher un dégradé blanc->rouge basé sur `Channel_Target_Display`.
        try:
            import numpy as np
            from matplotlib.colors import LinearSegmentedColormap

            if "Channel_Target_Display" in mesh.array_names:
                # Affichage ciblé : enlever la couche fantôme translucide
                # et utiliser un dégradé blanc→blanc→rouge concentré près de 1.0.
                # Le rendu principal reste la couche Target (0→1).

                # Colormap : blanc jusqu'à 0.80, légère montée puis rouge vif à 1.0
                cmap = LinearSegmentedColormap.from_list(
                    'white_to_red_sharp',
                    [
                        (0.0, '#ffffff'),
                        (0.70, "#ffffff"),
                        (0.80, "#f75c5c"),
                        (0.95, "#FF0800"),
                        (1.00, "#000000"),
                    ]
                )

                sbar = {
                    'title': 'Target (mm)',
                    'color': 'black',
                    'n_labels': 5,
                    'vertical': True,
                    'title_font_size': 12,
                    'label_font_size': 10,
                    'position_x': 0.85,
                    'position_y': 0.12,
                }

                self._add_tracked_mesh(si, mesh, "channels_target",
                    scalars="Channel_Target_Display",
                    cmap=cmap,
                    clim=[0.0, 1.0],
                    opacity=1.0,
                    scalar_bar_args=sbar,
                    show_scalar_bar=True)

                # Annoter la valeur cible si disponible
                if 'channel_target' in mesh.field_data:
                    try:
                        tgt = float(mesh.field_data['channel_target'][0])
                        key = self._view_key(si)
                        self.plotter.add_text(
                            f"Target: {tgt:.1f} mm",
                            position='upper_right', font_size=9,
                            color='black', name=f"{key}_ch_target_info")
                        if key not in self.view_actors:
                            self.view_actors[key] = []
                        self.view_actors[key].append(f"{key}_ch_target_info")
                    except Exception:
                        pass
                return

            # Fallback: render simple channels (no border)
            if "Channel_Region" not in mesh.array_names:
                self._add_tracked_mesh(si, mesh, "channels_fallback",
                    color="lightgrey", opacity=0.3, show_scalar_bar=False)
                return

            region = mesh["Channel_Region"]
            wt = mesh["Wall_Thickness"]

            display = np.zeros(len(region), dtype=float)
            display[region == 0] = -1.0  # rest
            display[region == 1] = wt[region == 1]  # channels
            mesh["Channel_Display"] = display

            # Simple colormap: gris pour le reste, bleu->red sur [0,5]
            colors_combined = [(-1.0, '#e0e0e0'), (0.0, '#2166ac'), (2.5, '#d73027'), (5.0, '#2166ac')]
            min_val, max_val = -1.0, 5.0
            norm_colors = [((val - min_val) / (max_val - min_val), col) for val, col in colors_combined]
            combined_cmap = LinearSegmentedColormap.from_list('channel_combined', norm_colors, N=512)
            sbar = {'title': 'Channels (mm)', 'color': 'black', 'n_labels': 6, 'vertical': True, 'title_font_size': 10, 'label_font_size': 8}
            self._add_tracked_mesh(si, mesh, "channels", scalars="Channel_Display", cmap=combined_cmap, clim=[min_val, max_val], opacity=1.0, scalar_bar_args=sbar, show_scalar_bar=True)
        except Exception:
            pass
    
    def _render_laplacian(self, si, mesh, mdef, lap_threshold=0.10, min_area_mm2=0.0):
        """Laplacien de l'épaisseur — visualisation binaire pour détecter les isthmes.
        Utilise le champ Laplacian déjà calculé par data_manager (pas de recalcul).
          - Bleu  = creux (∇²T ≥  seuil) → zones d'amincissement local
          - Rouge = bosses (∇²T ≤ -seuil) → zones d'épaississement local
          - Gris  = neutre (|∇²T| < seuil)
        Filtre de surface minimale via composantes connexes PyVista.
        """
        self.plotter.subplot(*si)

        if 'Laplacian' not in mesh.array_names:
            self._add_tracked_mesh(si, mesh, "lap_fallback",
                color="lightgrey", opacity=1.0, show_scalar_bar=False)
            return

        try:
            from matplotlib.colors import ListedColormap
            import numpy as np

            lap = mesh["Laplacian"].astype(float)

            # Classification ternaire
            category = np.zeros(len(lap), dtype=np.int8)
            category[lap >=  lap_threshold] = 1   # creux → bleu
            category[lap <= -lap_threshold] = 2   # bosses → rouge

            # --- Filtre surface minimale par composantes connexes ---
            if min_area_mm2 > 0:
                try:
                    for cat_val in (1, 2):
                        cat_mask = (category == cat_val)
                        if not np.any(cat_mask):
                            continue

                        # Extraire le sous-mesh de cette catégorie
                        ids = np.where(cat_mask)[0]
                        sub = mesh.extract_points(ids, adjacent_cells=True)
                        if sub.n_points == 0:
                            continue

                        # Composantes connexes
                        labeled = sub.connectivity(largest=False)
                        if 'RegionId' not in labeled.point_data:
                            continue

                        region_ids = labeled.point_data['RegionId']
                        n_regions = int(region_ids.max()) + 1

                        # Pour chaque région, estimer la surface et décider si on garde
                        # Trouver les points originaux à zapper
                        for r in range(n_regions):
                            r_mask = (region_ids == r)
                            sub_r = labeled.extract_points(
                                np.where(r_mask)[0], adjacent_cells=True)
                            try:
                                area = sub_r.area
                            except Exception:
                                area = sub_r.n_points  # fallback (nb points)
                            if area < min_area_mm2:
                                # Retrouver les indices originaux via proximité
                                r_pts = sub_r.points
                                if len(r_pts) == 0:
                                    continue
                                from scipy.spatial import cKDTree
                                tree = cKDTree(mesh.points)
                                _, orig_ids = tree.query(r_pts, k=1)
                                category[orig_ids] = 0  # → neutre
                except Exception:
                    pass

            mesh = mesh.copy()
            mesh["Laplacian_Category"] = category.astype(float)

            lap_cmap = ListedColormap(['#d0d0d0', '#d73027', '#2166ac'])
            sbar = {
                'title': 'Laplacian',
                'color': 'black',
                'n_labels': 3,
                'vertical': True,
                'title_font_size': 10,
                'label_font_size': 8,
            }
            self._add_tracked_mesh(si, mesh, "laplacian",
                scalars="Laplacian_Category",
                cmap=lap_cmap,
                clim=[0, 2],
                opacity=1.0,
                scalar_bar_args=sbar,
                show_scalar_bar=True)
        except Exception:
            pass
    
    def _render_isthmus(self, si, mesh, mdef):
        """Vue 5 de LV_topology — IDENTIQUE :
        p.add_mesh(mesh_gradient, color="lightgrey", opacity=0.1)
        p.add_mesh(isthme, color="purple", opacity=0.8)
        p.add_mesh(bounds, color="yellow", point_size=8, render_points_as_spheres=True)
        """
        self.plotter.subplot(*si)
        
        # Ghost mesh
        self._add_tracked_mesh(si, mesh, "isthmus_ghost",
            color="lightgrey", opacity=0.1, show_scalar_bar=False)
        
        # Isthmus zones (violet)
        try:
            isthme = mesh.threshold(0.5, scalars="Isthmus_Probability")
            if isthme.n_points > 0:
                self._add_tracked_mesh(si, isthme, "isthmus_zones",
                    color="purple", opacity=0.8, show_scalar_bar=False)
        except:
            pass
        
        # Boundaries (jaune, sphères)
        try:
            bounds = mesh.threshold(0.5, scalars="Isthmus_Boundaries")
            if bounds.n_points > 0:
                self._add_tracked_mesh(si, bounds, "isthmus_bounds",
                    color="yellow", point_size=8,
                    render_points_as_spheres=True,
                    show_scalar_bar=False)
        except:
            pass
    
    def _render_channelness(self, si, mesh, mdef):
        """Rendu channelness Cedilnik — ventricule entier coloré, isthmes en surbrillance"""
        self.plotter.subplot(*si)
        
        if "Channelness" not in mesh.array_names:
            self._add_tracked_mesh(si, mesh, "channelness_fallback",
                color="lightgrey", opacity=1.0, show_scalar_bar=False)
            return
        
        try:
            # Colormap ventricule entier (pas de seuil, pas de ghost)
            sbar = {
                'title': 'Channelness',
                'color': 'black',
                'vertical': True,
                'title_font_size': 11,
                'label_font_size': 9,
                'n_labels': 5,
                'position_x': 0.85,
                'position_y': 0.15,
            }
            
            self._add_tracked_mesh(si, mesh, "channelness_full",
                scalars="Channelness",
                cmap='inferno',
                clim=[0, 1],
                opacity=1.0,
                smooth_shading=True,
                show_edges=False,
                scalar_bar_args=sbar,
                show_scalar_bar=True)
            
            # Annotations : paramètres et statistiques
            info_parts = []
            if 'p_threshold' in mesh.field_data:
                p_val = float(mesh.field_data['p_threshold'][0])
                info_parts.append(f"p = {p_val:.1f} mm")
            if 'mean_delay' in mesh.field_data:
                md = float(mesh.field_data['mean_delay'][0])
                info_parts.append(f"Delay: {md:.2f}×")
            if 'n_pacing' in mesh.field_data:
                np_val = int(mesh.field_data['n_pacing'][0])
                info_parts.append(f"Pacing: {np_val} sites")
            if 'channelness_max' in mesh.field_data:
                mx = float(mesh.field_data['channelness_max'][0])
                info_parts.append(f"Max: {mx:.3f}")
            
            if info_parts:
                key = self._view_key(si)
                self.plotter.add_text(
                    "\n".join(info_parts),
                    position='upper_right',
                    font_size=9,
                    color='white',
                    shadow=True,
                    name=f"{key}_channelness_info"
                )
                if key not in self.view_actors:
                    self.view_actors[key] = []
                self.view_actors[key].append(f"{key}_channelness_info")
        
        except Exception:
            pass
    
    def _render_activation(self, si, mesh, mdef):
        """Comme simu_conduct: Carte isochronale avec cmap gist_rainbow + isochrones noires"""
        self.plotter.subplot(*si)
        sbar = {'title': 'Activation Time (ms)', 'color': 'black', 'n_labels': 6,
                'vertical': True, 'title_font_size': 12, 'label_font_size': 10,
                'position_x': 0.85, 'position_y': 0.15}
        
        # Déterminer la clim par défaut (max du mesh ou 750 ms par défaut)
        try:
            act = mesh["Activation_Time"]
            finite = act[~np.isinf(act)]
            if len(finite) > 0:
                clim_max = np.max(finite)
            else:
                clim_max = 750.0
        except:
            clim_max = 750.0
        
        # Colormap gist_rainbow identique à create_isochronal_map() de simu_conduct.py
        self._add_tracked_mesh(si, mesh, "activation",
            scalars="Activation_Time", cmap="gist_rainbow",
            smooth_shading=True, show_edges=False,
            clim=[0.0, clim_max],
            scalar_bar_args=sbar, show_scalar_bar=True)
        
        # Isochrones tous les 10ms (lignes noires) — identique simu_conduct.py
        ISOCHRONE_INTERVAL = 10.0
        try:
            act = mesh["Activation_Time"]
            finite = act[~np.isinf(act)]
            if len(finite) == 0:
                return
            max_time = np.max(finite)
            isochrone_times = np.arange(ISOCHRONE_INTERVAL, max_time, ISOCHRONE_INTERVAL)
            
            for k, iso_time in enumerate(isochrone_times):
                try:
                    contour = mesh.contour([iso_time], scalars="Activation_Time")
                    if contour.n_points > 0:
                        self._add_tracked_mesh(si, contour, f"iso_{k}",
                            color="black", line_width=2, opacity=0.6,
                            render_lines_as_tubes=False,
                            show_scalar_bar=False)
                except:
                    pass
        except:
            pass
    
    def _render_scar_burden(self, si, mesh, mdef):
        """Classification épaisseur 4 zones (colormap catégorielle) :
          0 = Sain   (T > healthy_thresh)       → vert
          1 = Border (border_thresh < T ≤ healthy_thresh)  → jaune
          2 = Scar   (dense_thresh < T ≤ border_thresh)  → orange
          3 = Dense  (T ≤ dense_thresh)       → rouge
        + grade Utah (bas gauche)
        Thresholds are read from mesh field_data if available.
        """
        self.plotter.subplot(*si)

        if "Scar_Burden_Display" not in mesh.array_names:
            self._add_tracked_mesh(si, mesh, "scar_burden_fallback",
                color="lightgrey", opacity=1.0, show_scalar_bar=False)
            return

        try:
            from matplotlib.colors import ListedColormap

            # 4 couleurs catégorielles
            scar_cmap = ListedColormap(['#4caf50', '#ffeb3b', '#ff9800', '#c62828'])

            sbar = {
                'title': 'Zone (WT)',
                'color': 'black',
                'n_labels': 0,
                'vertical': True,
                'title_font_size': 10,
                'label_font_size': 8,
            }

            self._add_tracked_mesh(si, mesh, "scar_burden",
                scalars="Scar_Burden_Display",
                cmap=scar_cmap,
                clim=[0, 3],
                opacity=1.0,
                scalar_bar_args=sbar,
                show_scalar_bar=True)

            # Overlay texte bas-gauche
            lines = []
            # Read thresholds from field_data for display
            t_healthy = 5.0
            t_border = 4.0
            t_dense = 2.0
            if 'thresh_healthy' in mesh.field_data:
                t_healthy = float(mesh.field_data['thresh_healthy'][0])
            if 'thresh_border' in mesh.field_data:
                t_border = float(mesh.field_data['thresh_border'][0])
            if 'thresh_dense' in mesh.field_data:
                t_dense = float(mesh.field_data['thresh_dense'][0])
            if 'wt_healthy_pct' in mesh.field_data:
                h = float(mesh.field_data['wt_healthy_pct'][0])
                b = float(mesh.field_data['wt_border_pct'][0])
                s = float(mesh.field_data['wt_scar_pct'][0])
                d = float(mesh.field_data['wt_dense_pct'][0])
                lines += [
                    f"Sain    >{t_healthy:.0f}mm   {h:.1f}%",
                    f"Border  {t_border:.0f}-{t_healthy:.0f}mm  {b:.1f}%",
                    f"Scar    {t_dense:.0f}-{t_border:.0f}mm  {s:.1f}%",
                    f"Dense   <{t_dense:.0f}mm   {d:.1f}%",
                ]
            if 'scar_burden_pct' in mesh.field_data:
                burden = float(mesh.field_data['scar_burden_pct'][0])
                if burden < 5:    ul = 'I - Minimal'
                elif burden < 20: ul = 'II - Mild'
                elif burden < 35: ul = 'III - Moderate'
                else:             ul = 'IV - Severe'
                lines.append(f"Utah {ul}")

            if lines:
                key = self._view_key(si)
                self.plotter.add_text(
                    "\n".join(lines),
                    position='upper_right',
                    font_size=9,
                    color='black',
                    shadow=False,
                    name=f"{key}_wt_info"
                )
                if key not in self.view_actors:
                    self.view_actors[key] = []
                self.view_actors[key].append(f"{key}_wt_info")

        except Exception:
            pass
    
    def _render_scar_distribution(self, si, mesh, mdef):
        """Comme advanced_scar_analysis: Scar segmentation by location"""
        self.plotter.subplot(*si)
        # mesh is a dict-like with endo/intra/epi — handled by caller
        self._add_tracked_mesh(si, mesh, "scar_dist",
            color='purple', opacity=1.0, show_scalar_bar=False)

    def _render_combined_score(self, si, mesh, mdef):
        """Rendu du score composite configurable [0, 1] avec plasma cmap
        + overlay des métriques impliquées et leurs poids."""
        self.plotter.subplot(*si)

        if 'Combined_Custom' not in mesh.array_names:
            self._add_tracked_mesh(si, mesh, 'combined_fallback',
                color='lightgrey', opacity=1.0, show_scalar_bar=False)
            return

        # Colormap discrète qualitative avec intervalles
        from matplotlib.colors import ListedColormap
        colors = [
            '#F5E7C6',  # Bleu clair (lisse)
            "#ffc76d",  # Orange (modéré)
            '#FF6D1F',  # Rouge (rugueux)
            "#dd2e1b",  # Violet (très rugueux)
            "#222222"   # Blanc (extrême)
        ]
        discrete_cmap = ListedColormap(colors)
        
        sbar = {
            'title': 'Score [0-1]',
            'color': 'black',
            'n_labels': 5,
            'vertical': True,
            'title_font_size': 10,
            'label_font_size': 8,
        }
        self._add_tracked_mesh(si, mesh, 'combined',
            scalars='Combined_Custom',
            cmap=discrete_cmap,
            clim=[0, 1],
            opacity=1.0,
            scalar_bar_args=sbar,
            show_scalar_bar=True)

        # Overlay : liste des composantes et leurs poids
        lines = ['SCORE COMBINÉ']
        if 'combined_config' in mesh.field_data:
            for entry in mesh.field_data['combined_config']:
                parts = str(entry).split('*')
                if len(parts) == 2:
                    try:
                        w = float(parts[1])
                        lines.append(f'  {parts[0]:<18s} {w*100:.0f}%')
                    except ValueError:
                        lines.append(f'  {entry}')
        if len(lines) > 1:
            key = self._view_key(si)
            self.plotter.add_text(
                '\n'.join(lines),
                position='upper_right',
                font_size=8,
                color='black',
                shadow=False,
                name=f'{key}_combined_info'
            )
            if key not in self.view_actors:
                self.view_actors[key] = []
            self.view_actors[key].append(f'{key}_combined_info')

    # -------------------------------------------------------------------------
    # ISOCHRONES OVERLAY + PACING SITE (appelé depuis main_window sur demande)
    # -------------------------------------------------------------------------
    ISOCHRONE_INTERVAL = 10.0  # ms

    def _resolve_activation(self, subplot_index: tuple, pacing_side: str, data_manager=None):
        """Résout le tableau d'activation et le mesh de référence pour RV ou LV.
        Retourne (act_arr, ref_mesh, ref_pts) ou (None, None, None)."""
        key = self._view_key(subplot_index)
        current_mesh = self.current_meshes.get(key)

        if pacing_side == 'RV':
            array_name = 'Activation_RV_Time'
            dm_keys = ('activation_RV', 'activation_rv', 'activation_RV')
        else:
            array_name = 'Activation_LV_Time'
            dm_keys = ('activation_LV', 'activation_lv', 'activation_LV')

        act_arr = None
        ref_mesh = None

        if current_mesh is not None and array_name in current_mesh.array_names:
            act_arr = current_mesh[array_name].astype(float)
            ref_mesh = current_mesh
        elif data_manager is not None:
            for dk in dm_keys:
                if dk in data_manager.computed_metrics:
                    src = data_manager.computed_metrics[dk]
                    act_arr = src['Activation_Time'].astype(float)
                    ref_mesh = src
                    break

        if act_arr is None or ref_mesh is None:
            return None, None, None

        if current_mesh is not None and len(act_arr) != current_mesh.n_points:
            tree = cKDTree(ref_mesh.points)
            _, idx = tree.query(current_mesh.points, k=1)
            act_arr = act_arr[idx]
            ref_mesh = current_mesh

        return act_arr, ref_mesh, ref_mesh.points

    def draw_isochrones(self, subplot_index: tuple, data_manager=None,
                        rv: bool = True, lv: bool = True):
        """Trace les isochrones par-dessus le subplot indiqué.
        rv / lv : activer chaque côté indépendamment."""
        i, j = subplot_index
        key = self._view_key(subplot_index)
        self.plotter.subplot(i, j)

        sides = []
        if rv:
            sides.append(('RV', '#1565c0'))
        if lv:
            sides.append(('LV', '#b71c1c'))

        for pacing_side, iso_color in sides:
            act_arr, ref_mesh, _ = self._resolve_activation(subplot_index, pacing_side, data_manager)
            if act_arr is None:
                continue

            finite = act_arr[~np.isinf(act_arr)]
            if len(finite) == 0:
                continue

            max_time = float(np.max(finite))
            iso_times = np.arange(self.ISOCHRONE_INTERVAL, max_time, self.ISOCHRONE_INTERVAL)
            tag = pacing_side.lower()

            work_mesh = ref_mesh.copy()
            work_mesh['_iso_time'] = act_arr

            if key not in self.view_actors:
                self.view_actors[key] = []

            for k_iso, iso_t in enumerate(iso_times):
                try:
                    contour = work_mesh.contour([iso_t], scalars='_iso_time')
                    if contour.n_points > 0:
                        iso_name = f'{key}_iso_{tag}_{k_iso}'
                        self.view_actors[key].append(iso_name)
                        self.plotter.add_mesh(
                            contour,
                            name=iso_name,
                            color=iso_color,
                            line_width=1.5,
                            opacity=0.75,
                            render_lines_as_tubes=False,
                            show_scalar_bar=False,
                        )
                except Exception:
                    pass

    def remove_isochrones(self, subplot_index: tuple, pacing_side: str = 'both'):
        """Supprime les acteurs isochrones du subplot.
        pacing_side : 'RV', 'LV' ou 'both'."""
        key = self._view_key(subplot_index)
        if key not in self.view_actors:
            return
        tags_to_remove = set()
        if pacing_side in ('RV', 'both'):
            tags_to_remove.add('_iso_rv_')
        if pacing_side in ('LV', 'both'):
            tags_to_remove.add('_iso_lv_')

        remaining = []
        for actor_name in self.view_actors[key]:
            remove = any(t in actor_name for t in tags_to_remove)
            if remove:
                try:
                    self.plotter.remove_actor(actor_name)
                except Exception:
                    pass
            else:
                remaining.append(actor_name)
        self.view_actors[key] = remaining

    # ── Pacing site marker ──────────────────────────────────────────────────

    def draw_pacing_site(self, subplot_index: tuple, pacing_side: str,
                         data_manager=None):
        """Affiche un triangle jaune au site de pacing (RV ou LV).
        Le point est lu depuis data_manager._pacing_points[pacing_side].
        Si la simulation n'a pas encore été lancée, on la lance."""
        i, j = subplot_index
        key = self._view_key(subplot_index)
        self.plotter.subplot(i, j)

        # Lancer la simulation si nécessaire pour obtenir le point
        if data_manager is not None:
            pacing_points = getattr(data_manager, '_pacing_points', {})
            if pacing_side not in pacing_points:
                try:
                    data_manager.simulate_activation(pacing_side)
                    pacing_points = getattr(data_manager, '_pacing_points', {})
                except Exception:
                    pass

            pt = pacing_points.get(pacing_side)
            if pt is None:
                return

            import pyvista as _pv
            pt_arr = np.array([pt])
            marker = _pv.PolyData(pt_arr)

            # Triangle glyphe (flèche vers le bas pour indiquer le site)
            cone = _pv.Cone(direction=(0, 0, -1), height=4.0, radius=2.5,
                            resolution=3)  # 3 faces → triangle

            marker_name = f'{key}_pacing_{pacing_side.lower()}'
            if key not in self.view_actors:
                self.view_actors[key] = []
            self.view_actors[key].append(marker_name)

            self.plotter.add_mesh(
                marker.glyph(geom=cone, scale=False),
                name=marker_name,
                color='yellow',
                opacity=1.0,
                show_scalar_bar=False,
            )

    def remove_pacing_site(self, subplot_index: tuple, pacing_side: str):
        """Supprime le marqueur de site de pacing pour ce subplot."""
        key = self._view_key(subplot_index)
        if key not in self.view_actors:
            return
        marker_name = f'{key}_pacing_{pacing_side.lower()}'
        remaining = []
        for actor_name in self.view_actors[key]:
            if actor_name == marker_name:
                try:
                    self.plotter.remove_actor(actor_name)
                except Exception:
                    pass
            else:
                remaining.append(actor_name)
        self.view_actors[key] = remaining

    def _render_generic(self, si, mesh, mdef):
        """Rendu générique pour métriques non spécialisées"""
        self.plotter.subplot(*si)
        sbar = {'title': f'{mdef.unit}' if mdef.unit else mdef.name,
                'color': 'black', 'n_labels': 4, 'vertical': True,
                'title_font_size': 10, 'label_font_size': 8}
        
        if mdef.scalar_name in mesh.array_names:
            self._add_tracked_mesh(si, mesh, "generic",
                scalars=mdef.scalar_name, cmap=mdef.cmap, clim=mdef.clim,
                smooth_shading=True, scalar_bar_args=sbar, show_scalar_bar=True)
        else:
            self._add_tracked_mesh(si, mesh, "generic",
                color='lightgray', opacity=0.5, show_scalar_bar=False)

    def render_empty(self, subplot_index: tuple, message: str = "No data"):
        """Rend une vue vide avec un message"""
        i, j = subplot_index
        key = self._view_key(subplot_index)
        
        self._clear_subplot(subplot_index)
        self.plotter.subplot(i, j)
        
        self.plotter.add_text(
            message,
            position='upper_left',
            font_size=12,
            color='gray',
            name=f"{key}_empty_text"
        )
    
    def set_view(self, subplot_index: tuple, view_type: str = 'anterior'):
        """Configure la vue d'un subplot.
        
        Le mesh VTK est déjà correctement orienté dans les coordonnées
        médicales du scanner. On positionne la caméra en conséquence.
        
        Conventions anatomiques :
          X+ = droite patient (Right)
          Y+ = supérieur (Superior)
          Z+ = antérieur (Anterior)
        """
        i, j = subplot_index
        self.plotter.subplot(i, j)
        
        key = self._view_key(subplot_index)
        mesh = self.current_meshes.get(key)
        
        if view_type == 'anterior':
            # Vue antérieure : caméra devant le patient (Z+), regard vers Z-
            # View-up = Y+ (supérieur)
            if mesh is not None:
                center = mesh.center
                # Distance = diagonale de la bounding box
                bounds = mesh.bounds
                diag = ((bounds[1]-bounds[0])**2 + (bounds[3]-bounds[2])**2 + (bounds[5]-bounds[4])**2) ** 0.5
                cam_pos = (center[0], center[1], center[2] + diag * 1.5)
                self.plotter.camera_position = [
                    cam_pos,              # position caméra
                    center,               # focal point
                    (0, 1, 0)             # view up (supérieur)
                ]
            else:
                self.plotter.view_isometric()
        elif view_type == 'posterior':
            if mesh is not None:
                center = mesh.center
                bounds = mesh.bounds
                diag = ((bounds[1]-bounds[0])**2 + (bounds[3]-bounds[2])**2 + (bounds[5]-bounds[4])**2) ** 0.5
                cam_pos = (center[0], center[1], center[2] - diag * 1.5)
                self.plotter.camera_position = [
                    cam_pos, center, (0, 1, 0)
                ]
            else:
                self.plotter.view_isometric()
        elif view_type == 'left':
            if mesh is not None:
                center = mesh.center
                bounds = mesh.bounds
                diag = ((bounds[1]-bounds[0])**2 + (bounds[3]-bounds[2])**2 + (bounds[5]-bounds[4])**2) ** 0.5
                cam_pos = (center[0] - diag * 1.5, center[1], center[2])
                self.plotter.camera_position = [
                    cam_pos, center, (0, 1, 0)
                ]
            else:
                self.plotter.view_isometric()
        elif view_type == 'right':
            if mesh is not None:
                center = mesh.center
                bounds = mesh.bounds
                diag = ((bounds[1]-bounds[0])**2 + (bounds[3]-bounds[2])**2 + (bounds[5]-bounds[4])**2) ** 0.5
                cam_pos = (center[0] + diag * 1.5, center[1], center[2])
                self.plotter.camera_position = [
                    cam_pos, center, (0, 1, 0)
                ]
            else:
                self.plotter.view_isometric()
        elif view_type == 'superior':
            if mesh is not None:
                center = mesh.center
                bounds = mesh.bounds
                diag = ((bounds[1]-bounds[0])**2 + (bounds[3]-bounds[2])**2 + (bounds[5]-bounds[4])**2) ** 0.5
                cam_pos = (center[0], center[1] + diag * 1.5, center[2])
                self.plotter.camera_position = [
                    cam_pos, center, (0, 0, -1)
                ]
            else:
                self.plotter.view_isometric()
        elif view_type == 'inferior':
            if mesh is not None:
                center = mesh.center
                bounds = mesh.bounds
                diag = ((bounds[1]-bounds[0])**2 + (bounds[3]-bounds[2])**2 + (bounds[5]-bounds[4])**2) ** 0.5
                cam_pos = (center[0], center[1] - diag * 1.5, center[2])
                self.plotter.camera_position = [
                    cam_pos, center, (0, 0, 1)
                ]
            else:
                self.plotter.view_isometric()
        elif view_type == 'isometric':
            self.plotter.view_isometric()
        elif view_type == 'xy':
            self.plotter.view_xy()
        elif view_type == 'xz':
            self.plotter.view_xz()
        elif view_type == 'yz':
            self.plotter.view_yz()
        elif view_type in ('lao', 'rao', 'll', 'rl'):
            # Vues obliques à 40°
            # Conventions : X+=right, Y+=superior, Z+=anterior
            import math
            ang = math.radians(40)
            c40, s40 = math.cos(ang), math.sin(ang)
            if mesh is not None:
                center = mesh.center
                bounds = mesh.bounds
                diag = ((bounds[1]-bounds[0])**2 + (bounds[3]-bounds[2])**2 + (bounds[5]-bounds[4])**2) ** 0.5
                d = diag * 1.5
                cx, cy, cz = center
                if view_type == 'lao':
                    # Left Anterior Oblique : tourne vers la gauche (X-) depuis AP
                    cam_pos = (cx - d*s40, cy, cz + d*c40)
                elif view_type == 'rao':
                    # Right Anterior Oblique : tourne vers la droite (X+) depuis AP
                    cam_pos = (cx + d*s40, cy, cz + d*c40)
                elif view_type == 'll':
                    # Left Lateral : vue depuis la gauche pure
                    cam_pos = (cx - d, cy, cz)
                else:  # rl
                    # Right Lateral : vue depuis la droite pure
                    cam_pos = (cx + d, cy, cz)
                self.plotter.camera_position = [cam_pos, center, (0, 1, 0)]
            else:
                self.plotter.view_isometric()
    
    def link_all_views(self):
        """Synchronise toutes les vues (rotation, zoom)"""
        self.plotter.link_views()

    def _render_anatomical_channelness(self, si, mesh, mdef):
        """Channelness anatomique géométrique — couloirs viables étroits.

        Rendu flat en un seul mesh (pas de superposition → pas de z-fighting).
        L'ensemble du ventricule est coloré par Anatomical_Channelness :
          ~0   = scar non-viable ou tissu large  (plasma sombre)
          >0   = corridor viable étroit          (plasma clair)
        """
        self.plotter.subplot(*si)

        if 'Anatomical_Channelness' not in mesh.array_names:
            self._add_tracked_mesh(si, mesh, "anat_ch_fallback",
                color="lightgrey", opacity=1.0, show_scalar_bar=False)
            return

        try:
            sbar = {
                'title': 'Anatomical Channelness',
                'color': 'black',
                'vertical': True,
                'title_font_size': 11,
                'label_font_size': 9,
                'n_labels': 5,
                'position_x': 0.85,
                'position_y': 0.15,
            }

            # Ventricule entier — un seul mesh, pas de sous-maillages superposés
            self._add_tracked_mesh(si, mesh, "anat_ch_full",
                scalars='Anatomical_Channelness',
                cmap='plasma',
                clim=[0, 1],
                opacity=1.0,
                smooth_shading=True,
                show_edges=False,
                scalar_bar_args=sbar,
                show_scalar_bar=True)

            # Squelette (axe médian) en supra-surfacique cyan
            if 'Is_Skeleton' in mesh.array_names:
                try:
                    skel = mesh.threshold(value=0.5, scalars='Is_Skeleton')
                    if skel.n_points > 0:
                        self._add_tracked_mesh(si, skel, "anat_ch_skeleton",
                            color="#00e5ff",
                            point_size=6,
                            style='points',
                            render_points_as_spheres=True,
                            opacity=1.0,
                            show_scalar_bar=False)
                except Exception:
                    pass

            # Annotations
            info = []
            fd = mesh.field_data
            if 'h_min'          in fd: info.append(f"h_min: {float(fd['h_min'][0]):.1f} mm")
            if 'max_width'      in fd: info.append(f"max W: {float(fd['max_width'][0]):.1f} mm")
            if 'narrow_pct'     in fd: info.append(f"Corridor: {float(fd['narrow_pct'][0]):.1f}%")
            if 'n_skeleton_pts' in fd: info.append(f"Skeleton: {int(fd['n_skeleton_pts'][0])} pts")
            if 'max_lw_raw'     in fd: info.append(f"L/W max: {float(fd['max_lw_raw'][0]):.1f}")
            if info:
                key = self._view_key(si)
                self.plotter.add_text("\n".join(info), position='upper_right',
                    font_size=9, color='white', shadow=True,
                    name=f"{key}_anat_ch_info")
                if key not in self.view_actors:
                    self.view_actors[key] = []
                self.view_actors[key].append(f"{key}_anat_ch_info")

        except Exception:
            pass

    def _render_cv_map(self, si, mesh, mdef):
        """Rendu carte de conduction velocity — sigmoïde de Cedilnik"""
        self.plotter.subplot(*si)
        
        if "CV_ms" not in mesh.array_names:
            self._add_tracked_mesh(si, mesh, "cv_fallback",
                color="lightgrey", opacity=1.0, show_scalar_bar=False)
            return
        
        try:
            sbar = {
                'title': 'CV (m/s)',
                'color': 'black',
                'vertical': True,
                'title_font_size': 11,
                'label_font_size': 9,
                'n_labels': 6,
                'position_x': 0.85,
                'position_y': 0.15,
            }
            
            # Ventricule entier coloré par CV
            self._add_tracked_mesh(si, mesh, "cv_full",
                scalars="CV_ms",
                cmap='RdYlGn',
                clim=[0, 0.6],
                opacity=1.0,
                smooth_shading=True,
                show_edges=False,
                scalar_bar_args=sbar,
                show_scalar_bar=True)
            
            # Contour à v = 0.1 m/s (frontière blocked/slow) et v = 0.3 (slow/normal)
            try:
                cv_vals = mesh["CV_ms"]
                for boundary_val, label, col in [
                    (0.05, "Block", "black"),
                    (0.18, "Slow", "darkred")
                ]:
                    if np.min(cv_vals) < boundary_val < np.max(cv_vals):
                        contour = mesh.contour([boundary_val], scalars="CV_ms")
                        if contour.n_points > 0:
                            self._add_tracked_mesh(si, contour, f"cv_contour_{label}",
                                color=col, line_width=2, opacity=0.7,
                                show_scalar_bar=False)
            except Exception:
                pass
            
            # Annotations
            key = self._view_key(si)
            cv_vals = mesh["CV_ms"]
            pct_blocked = float(np.mean(cv_vals < 0.05) * 100)
            pct_slow = float(np.mean((cv_vals >= 0.05) & (cv_vals < 0.18)) * 100)
            info = f"Blocked: {pct_blocked:.1f}%\nSlow: {pct_slow:.1f}%"
            self.plotter.add_text(info, position='lower_left', font_size=9,
                color='white', shadow=True, name=f"{key}_cv_info")
            if key not in self.view_actors:
                self.view_actors[key] = []
            self.view_actors[key].append(f"{key}_cv_info")
        
        except Exception:
            pass
