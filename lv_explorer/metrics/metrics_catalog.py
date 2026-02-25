"""
Catalogue centralisé de toutes les métriques visualisables
"""

from dataclasses import dataclass
from typing import Callable, Optional, Dict
import numpy as np


@dataclass
class MetricDefinition:
    """Définition d'une métrique visualisable"""
    
    id: str
    name: str
    category: str
    compute_func: str  # Nom de la méthode dans DataManager
    scalar_name: str   # Nom du champ scalaire dans le mesh
    cmap: str
    clim: Optional[tuple] = None
    unit: str = ""
    description: str = ""
    requires_threshold: bool = False
    title_color: str = "black"


class MetricsCatalog:
    """Catalogue de toutes les métriques disponibles"""
    
    def __init__(self):
        self.metrics = self._build_catalog()
        self.categories = self._build_categories()
    
    def _build_catalog(self) -> Dict[str, MetricDefinition]:
        """Construit le catalogue complet"""
        
        catalog = {}
        
        # === TOPOGRAPHIE ===
        catalog['thickness'] = MetricDefinition(
            id='thickness',
            name='Thickness Map',
            category='Topography',
            compute_func='compute_wall_thickness',
            scalar_name='EPI_Distance',
            cmap='jet_r',
            clim=(0, 5),
            unit='mm',
            description='Carte EPI distance (identique LV_topology View 1)'
        )
        
        # Suppression de la carte d'épaisseur pariétale (endo+epi)
        
        catalog['ciaccio'] = MetricDefinition(
            id='ciaccio',
            name='Wavefront Gradient',
            category='Topography',
            compute_func='compute_ciaccio_ratio',
            scalar_name='Ciaccio_Ratio_Display',
            cmap='jet_r',
            clim=(0, 0.5),
            unit='',
            description='Gradient d\'épaisseur du front d\'onde'
        )
        
        catalog['local_entropy'] = MetricDefinition(
            id='local_entropy',
            name='Local Entropy',
            category='Topography',
            compute_func='compute_local_entropy',
            scalar_name='Local_Entropy',
            cmap='inferno',
            clim=(0, 1),
            unit='',
            description='Entropie locale — mesure d\'irrégularité morphologique (zones chaotiques)'
        )
        
        catalog['tri'] = MetricDefinition(
            id='tri',
            name='TRI (Terrain Ruggedness)',
            category='Topography',
            compute_func='compute_tri',
            scalar_name='TRI',
            cmap='terrain',
            clim=None,
            unit='mm',
            description='Terrain Ruggedness Index — rugosité locale basée sur les variations d\'épaisseur avec voisins'
        )
        
        catalog['deceleration'] = MetricDefinition(
            id='deceleration',
            name='Deceleration Zones',
            category='Topography',
            compute_func='compute_deceleration_zones',
            scalar_name='DZ_Mask',
            cmap='Reds',
            clim=(0, 1),
            unit='',
            description='Zones de décélération (ρ > seuil)',
            requires_threshold=True,
            title_color='red'
        )
        
        catalog['channels'] = MetricDefinition(
            id='channels',
            name='Channels + Border Zone',
            category='Topography',
            compute_func='compute_channels',
            scalar_name='Channel_Score',
            cmap='viridis',
            clim=(0, 5),
            unit='mm',
            description='Canaux de conduction (0-5mm) + zone bordante ajustable',
            title_color='blue'
        )
        
        catalog['laplacian'] = MetricDefinition(
            id='laplacian',
            name='Laplacian (Surface Curvature)',
            category='Topography',
            compute_func='compute_laplacian',
            scalar_name='Laplacian',
            cmap='coolwarm',
            clim=(-0.5, 0.5),
            unit='mm⁻¹',
            description='Laplacien binaire — Bleu = creux, Rouge = bosses (seuil ajustable, filtre surface)',
            title_color='darkgreen',
            requires_threshold=True
        )

        catalog['isthmus'] = MetricDefinition(
            id='isthmus',
            name='Predicted VT Isthmus',
            category='Topography',
            compute_func='compute_isthmus_prediction',
            scalar_name='Isthmus_Probability',
            cmap='plasma',
            clim=(0, 1),
            unit='',
            description='Probabilité d\'isthme de tachycardie ventriculaire',
            title_color='purple'
        )
        
        catalog['channelness'] = MetricDefinition(
            id='channelness',
            name='Channelness (Cedilnik)',
            category='Topography',
            compute_func='compute_channelness',
            scalar_name='Channelness',
            cmap='inferno',
            clim=(0, 1),
            unit='',
            description='Carte de channelness simulation Eikonal (Cedilnik, EP-Europace 2018)',
            title_color='darkorange'
        )

        catalog['anatomical_channelness'] = MetricDefinition(
            id='anatomical_channelness',
            name='Anatomical Channelness (L/W)',
            category='Topography',
            compute_func='compute_anatomical_channelness',
            scalar_name='Anatomical_Channelness',
            cmap='plasma',
            clim=(0, 1),
            unit='',
            description=(
                'Corridors viables étroit — score géométrique pur (pas d\'EP). '
                'C(x) = 0.6 × (1 − w/W_max) + 0.4 × L/W normalisé. '
                'W = largeur locale via transformée de distance; L = longueur du squelette.'
            ),
            title_color='teal'
        )
        
        catalog['cv_map'] = MetricDefinition(
            id='cv_map',
            name='Conduction Velocity (Cedilnik)',
            category='Simulation',
            compute_func='compute_cv_map',
            scalar_name='CV_ms',
            cmap='RdYlGn',
            clim=(0, 0.6),
            unit='m/s',
            description='Vitesse de conduction via fonction de transfert sigmoïde '
                        '(Cedilnik, EP-Europace 2018 — v(x) = sigmoid(WT))',
            title_color='green'
        )
        
        # === SCAR BURDEN ===
        catalog['scar_burden'] = MetricDefinition(
            id='scar_burden',
            name='WT Classification (Utah)',
            category='Scar',
            compute_func='compute_scar_burden',
            scalar_name='Scar_Burden_Display',
            cmap='RdYlGn_r',
            clim=(0, 3),
            unit='',
            description='Classification 4 zones : Sain >5mm | Border 4-5mm | Scar 2-4mm | Dense ≤2mm (+ Utah grade)',
            title_color='darkred'
        )
        
        catalog['scar_distribution'] = MetricDefinition(
            id='scar_distribution',
            name='Scar Distribution',
            category='Scar',
            compute_func='get_scar_distribution',
            scalar_name='',
            cmap='Set1',
            clim=None,
            unit='',
            description='Distribution endo/intra/epi de la cicatrice'
        )
        
        # === SCORE COMBINÉ ===
        catalog['combined_score'] = MetricDefinition(
            id='combined_score',
            name='⊕ Combined Score',
            category='Combined',
            compute_func='compute_custom_combined',
            scalar_name='Combined_Custom',
            cmap='plasma',
            clim=(0, 1),
            unit='',
            description='Score composite configurable : somme pondérée de métriques '
                        'normalisées (canal, gradient, channelness, entropie…)',
            title_color='darkviolet'
        )

        # === SIMULATION ===
        catalog['simulation'] = MetricDefinition(
            id='simulation',
            name='Simulation',
            category='Simulation',
            compute_func='simulate_activation',
            scalar_name='Activation_Time',
            cmap='turbo',
            clim=None,
            unit='ms',
            description='Carte d\'activation — clic droit pour choisir le site de pacing'
        )

        return catalog
    
    def _build_categories(self) -> Dict[str, list]:
        """Organise les métriques par catégorie"""
        
        categories = {}
        for metric_id, metric in self.metrics.items():
            if metric.category not in categories:
                categories[metric.category] = []
            categories[metric.category].append(metric_id)
        
        return categories
    
    def get_metric(self, metric_id: str) -> Optional[MetricDefinition]:
        """Récupère une métrique par son ID"""
        return self.metrics.get(metric_id)
    
    def get_by_category(self, category: str) -> list:
        """Récupère toutes les métriques d'une catégorie"""
        return self.categories.get(category, [])
    
    def get_all_ids(self) -> list:
        """Retourne la liste de tous les IDs"""
        return list(self.metrics.keys())
    
    def get_all_names(self) -> list:
        """Retourne la liste de tous les noms"""
        return [m.name for m in self.metrics.values()]