"""
Interface principale de LV Explorer
"""

import sys
import math
import pyvista as pv
from pyvistaqt import QtInteractor
from qtpy import QtWidgets, QtCore, QtGui
from qtpy.QtCore import QPropertyAnimation, QEasingCurve
import numpy as np

from ..core.data_manager import DataManager
from ..core.visualization_manager import VisualizationManager
from ..core.orientation_widget import HumanBustOrientationWidget
from ..metrics.metrics_catalog import MetricsCatalog
from .combined_score_dialog import CombinedScoreDialog

try:
    from ..core.dashboard_manager import DashboardWindow
    _HAS_DASHBOARD = True
except ImportError:
    _HAS_DASHBOARD = False


class _VizOverlayManager(QtCore.QObject):
    """Repositionne le widget overlay lors des redimensionnements du parent."""
    def __init__(self, overlay: QtWidgets.QWidget, target: QtWidgets.QWidget):
        super().__init__(target)
        self._overlay = overlay
        target.installEventFilter(self)

    def eventFilter(self, obj, event):
        if event.type() == QtCore.QEvent.Resize:
            self._reposition(obj)
        return False

    def _reposition(self, target: QtWidgets.QWidget):
        w = self._overlay.sizeHint().width()
        h = self._overlay.sizeHint().height()
        margin = 10
        x = margin
        y = max(margin, target.height() - h - margin)
        self._overlay.setGeometry(x, y, w, h)
        self._overlay.raise_()


def _optimal_grid(n: int) -> tuple:
    """
    Calcule la grille optimale (rows, cols) pour n vues.
    1→(1,1), 2→(1,2), 3→(1,3), 4→(2,2), 5→(2,3), 6→(2,3),
    7→(3,3), 8→(3,3), 9→(3,3), etc.
    """
    if n <= 0:
        return (1, 1)
    if n == 1:
        return (1, 1)
    if n == 2:
        return (1, 2)
    if n == 3:
        return (1, 3)
    if n == 4:
        return (2, 2)
    cols = math.ceil(math.sqrt(n))
    rows = math.ceil(n / cols)
    return (rows, cols)


class LVExplorerApp(QtWidgets.QMainWindow):
    """
    Application principale - Interface avec grille de visualisation adaptive
    """
    
    MAX_VIEWS = 9   # nombre maximum de vues simultanées
    
    def __init__(self, data_path: str = None):
        super().__init__()
        
        self.data_path = data_path
        self.data_manager = None
        self.viz_manager = None
        self.plotter = None
        self.orientation_widget = None
        self.catalog = MetricsCatalog()
        self._combined_configs = {}   # view_idx → config list pour combined_score
        
        # Nombre de vues actuel et shape courante
        self._n_views = 1
        self._grid_shape = (1, 1)
        self._initialized_views = set()   # subplots (i,j) ayant déjà un métrique
        self._panel_collapsed = False
        self._panel_anim = None   # garder une référence pour éviter le GC
        
        # Configuration de la fenêtre
        self.setWindowTitle("LV Explorer - Ventricular Analysis")
        self.setGeometry(100, 100, 1600, 900)
        
        self._setup_ui()
        
        if data_path:
            self.load_patient(data_path)
    
    def _setup_ui(self):
        """Configure l'interface utilisateur"""
        
        # Widget central
        central_widget = QtWidgets.QWidget()
        self.setCentralWidget(central_widget)
        
        # Layout principal horizontal
        main_layout = QtWidgets.QHBoxLayout(central_widget)
        main_layout.setContentsMargins(0, 0, 0, 0)
        main_layout.setSpacing(0)
        
        # === PANNEAU GAUCHE (Contrôles) ===
        self.left_panel = self._create_left_panel()
        main_layout.addWidget(self.left_panel, stretch=0)

        # === BOUTON DE REPLIAGE DU PANNEAU ===
        self._panel_toggle_btn = QtWidgets.QPushButton("◄")
        self._panel_toggle_btn.setFixedWidth(16)
        self._panel_toggle_btn.setSizePolicy(
            QtWidgets.QSizePolicy.Fixed, QtWidgets.QSizePolicy.Expanding)
        self._panel_toggle_btn.setToolTip("Masquer / Afficher le panneau")
        self._panel_toggle_btn.setStyleSheet("""
            QPushButton {
                background: #e8e8e8; border: none; border-left: 1px solid #d0d0d0;
                color: #666; font-size: 11px; border-radius: 0;
            }
            QPushButton:hover { background: #1976D2; color: white; }
        """)
        self._panel_toggle_btn.clicked.connect(self._toggle_left_panel)
        main_layout.addWidget(self._panel_toggle_btn, stretch=0)
        
        # === ZONE CENTRALE (Visualisation) ===
        self.viz_container = QtWidgets.QWidget()
        self.viz_layout = QtWidgets.QVBoxLayout(self.viz_container)
        self.viz_layout.setContentsMargins(0, 0, 0, 0)
        self.viz_layout.setSpacing(0)
        main_layout.addWidget(self.viz_container, stretch=1)

        self._create_plotter(self._grid_shape)
        self._create_view_buttons_bar()

        # === BARRE DE MENU ===
        self._create_menu_bar()
    
    def _create_left_panel(self) -> QtWidgets.QWidget:
        """Crée le panneau latéral gauche — design moderne, minimaliste, sobre"""
        
        panel = QtWidgets.QWidget()
        panel.setMaximumWidth(420)
        panel.setMinimumWidth(340)
        panel.setStyleSheet("""
            QWidget {
                background-color: #fafafa;
                font-family: 'Segoe UI', 'Helvetica Neue', Arial, sans-serif;
                font-size: 13px;
                color: #222;
            }
            QGroupBox {
                font-weight: 600;
                font-size: 14px;
                border: 1px solid #e0e0e0;
                border-radius: 8px;
                margin-top: 12px;
                padding: 14px 10px 10px 10px;
                background-color: #ffffff;
            }
            QGroupBox::title {
                subcontrol-origin: margin;
                left: 12px;
                padding: 0 6px;
                color: #333;
            }
            QPushButton {
                background-color: #1976D2;
                color: white;
                border: none;
                padding: 10px 18px;
                border-radius: 6px;
                font-size: 13px;
                font-weight: 600;
                min-height: 28px;
            }
            QPushButton:hover {
                background-color: #1565C0;
            }
            QPushButton:pressed {
                background-color: #0D47A1;
            }
            QComboBox {
                font-size: 13px;
                padding: 6px 10px;
                border: 1px solid #ccc;
                border-radius: 6px;
                background: white;
                min-height: 26px;
            }
            QComboBox:hover {
                border-color: #1976D2;
            }
            QComboBox::drop-down {
                border: none;
                width: 24px;
            }
            QLabel {
                font-size: 13px;
            }
            QDoubleSpinBox, QSpinBox {
                font-size: 13px;
                padding: 4px 8px;
                border: 1px solid #ccc;
                border-radius: 6px;
                background: white;
                min-height: 24px;
            }
            QSlider::groove:horizontal {
                height: 6px;
                background: #e0e0e0;
                border-radius: 3px;
            }
            QSlider::handle:horizontal {
                background: #1976D2;
                width: 18px;
                height: 18px;
                margin: -6px 0;
                border-radius: 9px;
            }
            QSlider::handle:horizontal:hover {
                background: #1565C0;
            }
            QTextEdit {
                font-family: 'Consolas', 'Courier New', monospace;
                font-size: 11px;
                background-color: #ffffff;
                border: 1px solid #e0e0e0;
                border-radius: 6px;
                padding: 6px;
            }
            QCheckBox {
                font-size: 13px;
                spacing: 8px;
            }
            QCheckBox::indicator {
                width: 18px;
                height: 18px;
            }
        """)
        
        layout = QtWidgets.QVBoxLayout(panel)
        layout.setContentsMargins(10, 10, 10, 10)
        layout.setSpacing(10)
        
        # === 1. PATIENT INFO ===
        patient_group = QtWidgets.QGroupBox("Patient Data")
        patient_layout = QtWidgets.QVBoxLayout()
        
        self.patient_label = QtWidgets.QLabel("No patient loaded")
        self.patient_label.setStyleSheet("font-size: 13px; color: #888;")
        patient_layout.addWidget(self.patient_label)
        
        load_btn = QtWidgets.QPushButton("Load Patient...")
        load_btn.clicked.connect(self._on_load_patient)
        patient_layout.addWidget(load_btn)
        
        patient_group.setLayout(patient_layout)
        layout.addWidget(patient_group)
        
        # === 2. GRID / VIEW COUNT ===
        grid_group = QtWidgets.QGroupBox("Grid Layout")
        grid_layout = QtWidgets.QHBoxLayout()
        
        minus_btn = QtWidgets.QPushButton("−")
        minus_btn.setFixedSize(40, 40)
        minus_btn.setStyleSheet("font-size: 20px; font-weight: bold; padding: 0;")
        minus_btn.setToolTip("Retirer une vue")
        minus_btn.clicked.connect(self._on_remove_view)
        grid_layout.addWidget(minus_btn)
        
        self.view_count_label = QtWidgets.QLabel(str(self._n_views))
        self.view_count_label.setAlignment(QtCore.Qt.AlignCenter)
        self.view_count_label.setStyleSheet("font-size: 18px; font-weight: bold; min-width: 30px;")
        grid_layout.addWidget(self.view_count_label)
        
        plus_btn = QtWidgets.QPushButton("+")
        plus_btn.setFixedSize(40, 40)
        plus_btn.setStyleSheet("font-size: 20px; font-weight: bold; padding: 0;")
        plus_btn.setToolTip("Ajouter une vue")
        plus_btn.clicked.connect(self._on_add_view)
        grid_layout.addWidget(plus_btn)
        
        self.grid_label = QtWidgets.QLabel("(1×1)")
        self.grid_label.setStyleSheet("font-size: 12px; color: #888;")
        grid_layout.addWidget(self.grid_label)
        
        auto_btn = QtWidgets.QPushButton("Auto")
        auto_btn.setToolTip("Ajuste automatiquement le nombre de vues\nau nombre de métriques sélectionnées")
        auto_btn.setMaximumWidth(60)
        auto_btn.clicked.connect(self._auto_adjust_views)
        grid_layout.addWidget(auto_btn)
        
        grid_group.setLayout(grid_layout)
        layout.addWidget(grid_group)
        
        # === 3. VIEW CONFIGURATION (dynamic combos) ===
        self.views_group = QtWidgets.QGroupBox("View Configuration")
        self.views_inner_layout = QtWidgets.QVBoxLayout()
        self.views_group.setLayout(self.views_inner_layout)
        layout.addWidget(self.views_group)
        
        self.view_combos = {}
        self._rebuild_view_combos()
        
        # === 4. THICKNESS SLIDER (visible uniquement pour thickness / parietal_thickness) ===
        self.thickness_group = QtWidgets.QGroupBox("Épaisseur pariétale — Échelle")
        thickness_layout = QtWidgets.QVBoxLayout()

        self.thickness_label = QtWidgets.QLabel("Max : 6.0 mm  (6 zones × 1.00 mm)")
        self.thickness_label.setStyleSheet("font-size: 13px; font-weight: bold;")
        thickness_layout.addWidget(self.thickness_label)

        self.thickness_slider = QtWidgets.QSlider(QtCore.Qt.Horizontal)
        self.thickness_slider.setMinimum(10)   # 1.0 mm
        self.thickness_slider.setMaximum(150)  # 15.0 mm
        self.thickness_slider.setValue(60)     # 6.0 mm par défaut
        self.thickness_slider.setTickInterval(10)
        self.thickness_slider.setTickPosition(QtWidgets.QSlider.TicksBelow)
        self.thickness_slider.valueChanged.connect(self._on_thickness_slider_changed)
        thickness_layout.addWidget(self.thickness_slider)

        thickness_hint = QtWidgets.QLabel("6 couleurs discrètes — de 0 mm au max")
        thickness_hint.setStyleSheet("font-size: 11px; color: #888;")
        thickness_layout.addWidget(thickness_hint)

        self.thickness_group.setLayout(thickness_layout)
        self.thickness_group.setVisible(False)
        layout.addWidget(self.thickness_group)

        # === 4b. DZ THRESHOLD SLIDER (visible uniquement pour Deceleration Zones) ===
        self.dz_group = QtWidgets.QGroupBox("DZ Threshold")
        dz_layout = QtWidgets.QVBoxLayout()
        
        self.dz_label = QtWidgets.QLabel("Threshold: 0.33")
        self.dz_label.setStyleSheet("font-size: 13px; font-weight: bold;")
        dz_layout.addWidget(self.dz_label)
        
        self.dz_slider = QtWidgets.QSlider(QtCore.Qt.Horizontal)
        self.dz_slider.setMinimum(1)    # 0.01
        self.dz_slider.setMaximum(100)  # 1.00
        self.dz_slider.setValue(33)     # 0.33 par défaut
        self.dz_slider.setTickInterval(10)
        self.dz_slider.setTickPosition(QtWidgets.QSlider.TicksBelow)
        self.dz_slider.valueChanged.connect(self._on_dz_threshold_changed)
        dz_layout.addWidget(self.dz_slider)
        
        dz_hint = QtWidgets.QLabel("Ciaccio ratio ρ ≥ threshold → DZ")
        dz_hint.setStyleSheet("font-size: 11px; color: #888;")
        dz_layout.addWidget(dz_hint)
        
        self.dz_group.setLayout(dz_layout)
        self.dz_group.setVisible(False)  # Initialement caché
        layout.addWidget(self.dz_group)

        # === 4c. BORDER RADIUS SLIDER (visible uniquement pour Channels) ===
        self.border_group = QtWidgets.QGroupBox("Border Zone Radius")
        border_layout = QtWidgets.QVBoxLayout()
        
        self.border_label = QtWidgets.QLabel("Radius: 10.0 mm")
        self.border_label.setStyleSheet("font-size: 13px; font-weight: bold;")
        border_layout.addWidget(self.border_label)
        
        self.border_slider = QtWidgets.QSlider(QtCore.Qt.Horizontal)
        self.border_slider.setMinimum(0)     # 0 mm
        self.border_slider.setMaximum(50)    # 50 mm
        self.border_slider.setValue(10)      # 10 mm par défaut
        self.border_slider.setTickInterval(5)
        self.border_slider.setTickPosition(QtWidgets.QSlider.TicksBelow)
        self.border_slider.valueChanged.connect(self._on_border_radius_changed)
        border_layout.addWidget(self.border_slider)
        
        border_hint = QtWidgets.QLabel("Distance autour des canaux (0-5mm zone)")
        border_hint.setStyleSheet("font-size: 11px; color: #888;")
        border_layout.addWidget(border_hint)
        
        self.border_group.setLayout(border_layout)
        self.border_group.setVisible(False)
        layout.addWidget(self.border_group)

        # === 4d. LAPLACIAN THRESHOLD SLIDER + SURFACE FILTER ===
        self.lap_group = QtWidgets.QGroupBox("Laplacian")
        lap_layout = QtWidgets.QVBoxLayout()

        # --- Seuil ∇²T ---
        lap_thresh_lbl = QtWidgets.QLabel("Seuil |∇²T| :")
        lap_thresh_lbl.setStyleSheet("font-size: 13px; font-weight: bold;")
        lap_layout.addWidget(lap_thresh_lbl)

        self.lap_thresh_spin = QtWidgets.QDoubleSpinBox()
        self.lap_thresh_spin.setRange(0.001, 200.0)
        self.lap_thresh_spin.setValue(0.10)
        self.lap_thresh_spin.setSingleStep(0.01)
        self.lap_thresh_spin.setDecimals(3)
        self.lap_thresh_spin.setKeyboardTracking(False)
        self.lap_thresh_spin.valueChanged.connect(self._on_lap_threshold_changed)
        lap_layout.addWidget(self.lap_thresh_spin)

        lap_hint = QtWidgets.QLabel("Bleu = creux (∇²T ≥ seuil)\nRouge = bosses (∇²T ≤ −seuil)")
        lap_hint.setStyleSheet("font-size: 11px; color: #888;")
        lap_layout.addWidget(lap_hint)

        # --- Surface minimale ---
        lap_area_lbl = QtWidgets.QLabel("Surface min (mm²) :")
        lap_area_lbl.setStyleSheet("font-size: 13px; font-weight: bold;")
        lap_layout.addWidget(lap_area_lbl)

        self.lap_area_spin = QtWidgets.QDoubleSpinBox()
        self.lap_area_spin.setRange(0.0, 5000.0)
        self.lap_area_spin.setValue(0.0)
        self.lap_area_spin.setSingleStep(10.0)
        self.lap_area_spin.setDecimals(0)
        self.lap_area_spin.setSuffix(" mm²")
        self.lap_area_spin.setKeyboardTracking(False)
        self.lap_area_spin.valueChanged.connect(self._on_lap_threshold_changed)
        lap_layout.addWidget(self.lap_area_spin)

        lap_area_hint = QtWidgets.QLabel("Ne garder que les zones > surface min")
        lap_area_hint.setStyleSheet("font-size: 11px; color: #888;")
        lap_layout.addWidget(lap_area_hint)

        self.lap_group.setLayout(lap_layout)
        self.lap_group.setVisible(False)
        layout.addWidget(self.lap_group)

        # === 4e. CEDILNIK DILATION RADIUS SLIDER ===
        self.cedilnik_group = QtWidgets.QGroupBox("Cedilnik BZ Dilation")
        cedilnik_layout = QtWidgets.QVBoxLayout()
        
        self.cedilnik_label = QtWidgets.QLabel("Radius: 2.0 mm")
        self.cedilnik_label.setStyleSheet("font-size: 13px; font-weight: bold;")
        cedilnik_layout.addWidget(self.cedilnik_label)
        
        self.cedilnik_slider = QtWidgets.QSlider(QtCore.Qt.Horizontal)
        self.cedilnik_slider.setMinimum(1)     # 0.5 mm (×2)
        self.cedilnik_slider.setMaximum(20)    # 10.0 mm (×2)
        self.cedilnik_slider.setValue(4)       # 2.0 mm par défaut (×2)
        self.cedilnik_slider.setTickInterval(2)
        self.cedilnik_slider.setTickPosition(QtWidgets.QSlider.TicksBelow)
        self.cedilnik_slider.valueChanged.connect(self._on_cedilnik_radius_changed)
        cedilnik_layout.addWidget(self.cedilnik_slider)
        
        cedilnik_hint = QtWidgets.QLabel("Intersection scar⊕R ∩ healthy⊕R (Cedilnik)")
        cedilnik_hint.setStyleSheet("font-size: 11px; color: #888;")
        cedilnik_layout.addWidget(cedilnik_hint)
        
        self.cedilnik_group.setLayout(cedilnik_layout)
        self.cedilnik_group.setVisible(False)
        layout.addWidget(self.cedilnik_group)

        # === 4f. CHANNELNESS p SLIDER (Cedilnik) ===
        self.channelness_group = QtWidgets.QGroupBox("Channelness (Cedilnik)")
        channelness_layout = QtWidgets.QVBoxLayout()
        
        self.channelness_label = QtWidgets.QLabel("Seuil p : 3.0 mm")
        self.channelness_label.setStyleSheet("font-size: 13px; font-weight: bold;")
        channelness_layout.addWidget(self.channelness_label)
        
        self.channelness_slider = QtWidgets.QSlider(QtCore.Qt.Horizontal)
        self.channelness_slider.setMinimum(15)    # 1.5 mm (×10)
        self.channelness_slider.setMaximum(50)    # 5.0 mm (×10)
        self.channelness_slider.setValue(30)      # 3.0 mm par défaut (×10)
        self.channelness_slider.setTickInterval(5)
        self.channelness_slider.setSingleStep(5)  # pas de 0.5 mm
        self.channelness_slider.setPageStep(5)
        self.channelness_slider.setTickPosition(QtWidgets.QSlider.TicksBelow)
        self.channelness_slider.valueChanged.connect(self._on_channelness_params_changed)
        channelness_layout.addWidget(self.channelness_slider)

        channelness_hint = QtWidgets.QLabel("Point d'inflexion sigmoïde T→CV (Cedilnik)")
        channelness_hint.setStyleSheet("font-size: 11px; color: #888;")
        channelness_layout.addWidget(channelness_hint)

        # ── Slider r : raideur du ralentissement exponentiel ──────────
        self.channelness_r_label = QtWidgets.QLabel("Raideur r : 2.0 mm⁻¹")
        self.channelness_r_label.setStyleSheet("font-size: 13px; font-weight: bold;")
        channelness_layout.addWidget(self.channelness_r_label)

        self.channelness_r_slider = QtWidgets.QSlider(QtCore.Qt.Horizontal)
        self.channelness_r_slider.setMinimum(5)   # 0.5 mm⁻¹
        self.channelness_r_slider.setMaximum(40)  # 4.0 mm⁻¹
        self.channelness_r_slider.setValue(20)    # 2.0 mm⁻¹ par défaut
        self.channelness_r_slider.setTickInterval(5)
        self.channelness_r_slider.setSingleStep(5)  # pas de 0.5 mm⁻¹
        self.channelness_r_slider.setPageStep(5)
        self.channelness_r_slider.setTickPosition(QtWidgets.QSlider.TicksBelow)
        self.channelness_r_slider.valueChanged.connect(self._on_channelness_params_changed)
        channelness_layout.addWidget(self.channelness_r_slider)

        r_hint = QtWidgets.QLabel("Raideur exposant. de la transition WT→CV")
        r_hint.setStyleSheet("font-size: 11px; color: #888;")
        channelness_layout.addWidget(r_hint)

        self.channelness_group.setLayout(channelness_layout)
        self.channelness_group.setVisible(False)
        layout.addWidget(self.channelness_group)

        # === 4g. ANATOMICAL CHANNELNESS (h_min + max_width) ===
        self.anat_ch_group = QtWidgets.QGroupBox("Anatomical Channelness")
        anat_ch_layout = QtWidgets.QVBoxLayout()

        # h_min slider
        anat_ch_layout.addWidget(QtWidgets.QLabel("Viable threshold (hₛᵢₙ)"))
        self.anat_ch_hmin_label = QtWidgets.QLabel("h_min: 3.0 mm")
        self.anat_ch_hmin_label.setObjectName("label-value")
        anat_ch_layout.addWidget(self.anat_ch_hmin_label)
        self.anat_ch_hmin_slider = QtWidgets.QSlider(QtCore.Qt.Horizontal)
        self.anat_ch_hmin_slider.setMinimum(10)   # 1.0 mm
        self.anat_ch_hmin_slider.setMaximum(70)   # 7.0 mm
        self.anat_ch_hmin_slider.setValue(30)     # 3.0 mm default
        self.anat_ch_hmin_slider.setTickInterval(10)
        self.anat_ch_hmin_slider.setTickPosition(QtWidgets.QSlider.TicksBelow)
        self.anat_ch_hmin_slider.valueChanged.connect(self._on_anat_ch_changed)
        anat_ch_layout.addWidget(self.anat_ch_hmin_slider)
        hmin_hint = QtWidgets.QLabel("T > h_min → viable tissue")
        hmin_hint.setObjectName("label-muted")
        anat_ch_layout.addWidget(hmin_hint)

        # max_width slider
        anat_ch_layout.addWidget(QtWidgets.QLabel("Max corridor width"))
        self.anat_ch_maxw_label = QtWidgets.QLabel("max W: 8.0 mm")
        self.anat_ch_maxw_label.setObjectName("label-value")
        anat_ch_layout.addWidget(self.anat_ch_maxw_label)
        self.anat_ch_maxw_slider = QtWidgets.QSlider(QtCore.Qt.Horizontal)
        self.anat_ch_maxw_slider.setMinimum(20)   # 2.0 mm
        self.anat_ch_maxw_slider.setMaximum(150)  # 15.0 mm
        self.anat_ch_maxw_slider.setValue(80)     # 8.0 mm default
        self.anat_ch_maxw_slider.setTickInterval(10)
        self.anat_ch_maxw_slider.setTickPosition(QtWidgets.QSlider.TicksBelow)
        self.anat_ch_maxw_slider.valueChanged.connect(self._on_anat_ch_changed)
        anat_ch_layout.addWidget(self.anat_ch_maxw_slider)
        maxw_hint = QtWidgets.QLabel("Corridors: viable AND width < max W")
        maxw_hint.setObjectName("label-muted")
        anat_ch_layout.addWidget(maxw_hint)

        self.anat_ch_group.setLayout(anat_ch_layout)
        self.anat_ch_group.setVisible(False)
        layout.addWidget(self.anat_ch_group)

        # === 4h. SIMULATION ACTIVATION (scar_decay) ===
        self.activation_group = QtWidgets.QGroupBox("Simulation d'activation")
        activation_layout = QtWidgets.QVBoxLayout()

        self.scar_decay_label = QtWidgets.QLabel("Decay cicatrice : 5.0 mm")
        self.scar_decay_label.setStyleSheet("font-size: 13px; font-weight: bold;")
        activation_layout.addWidget(self.scar_decay_label)

        self.scar_decay_slider = QtWidgets.QSlider(QtCore.Qt.Horizontal)
        self.scar_decay_slider.setMinimum(2)   # 1.0 mm  (×2)
        self.scar_decay_slider.setMaximum(40)  # 20.0 mm (×2)
        self.scar_decay_slider.setValue(10)    # 5.0 mm  (×2, défaut)
        self.scar_decay_slider.setSingleStep(1)   # 0.5 mm
        self.scar_decay_slider.setPageStep(2)
        self.scar_decay_slider.setTickInterval(2)
        self.scar_decay_slider.setTickPosition(QtWidgets.QSlider.TicksBelow)
        self.scar_decay_slider.valueChanged.connect(self._on_scar_decay_changed)
        activation_layout.addWidget(self.scar_decay_slider)

        decay_hint = QtWidgets.QLabel(
            "Portée du ralentissement exponentiel\n"
            "depuis la cicatrice (mm) — exp(−d/decay)")
        decay_hint.setStyleSheet("font-size: 11px; color: #888;")
        activation_layout.addWidget(decay_hint)

        # === Slider Activation Time Scale (clim pour la carte d'activation) ===
        self.at_scale_label = QtWidgets.QLabel("Échelle latence : 750 ms")
        self.at_scale_label.setStyleSheet("font-size: 13px; font-weight: bold;")
        activation_layout.addWidget(self.at_scale_label)

        self.at_scale_slider = QtWidgets.QSlider(QtCore.Qt.Horizontal)
        self.at_scale_slider.setMinimum(100)  # 100 ms
        self.at_scale_slider.setMaximum(1500) # 1500 ms
        self.at_scale_slider.setValue(750)    # 750 ms par défaut
        self.at_scale_slider.setSingleStep(10)
        self.at_scale_slider.setPageStep(50)
        self.at_scale_slider.setTickInterval(100)
        self.at_scale_slider.setTickPosition(QtWidgets.QSlider.TicksBelow)
        self.at_scale_slider.valueChanged.connect(self._on_at_scale_changed)
        activation_layout.addWidget(self.at_scale_slider)

        at_hint = QtWidgets.QLabel("Plage colorée de la carte : 0 à max (ms)")
        at_hint.setStyleSheet("font-size: 11px; color: #888;")
        activation_layout.addWidget(at_hint)

        # Pacing info label
        self.pacing_info_label = QtWidgets.QLabel("Clic droit sur le mesh pour\nchoisir le site de pacing")
        self.pacing_info_label.setStyleSheet("font-size: 12px; color: #1976D2; font-weight: 600;")
        activation_layout.addWidget(self.pacing_info_label)

        self.activation_group.setLayout(activation_layout)
        self.activation_group.setVisible(False)
        layout.addWidget(self.activation_group)

        # === 4i. CLASSIFICATION THRESHOLDS (Scar Burden) ===
        self.classif_group = QtWidgets.QGroupBox("Classification Thresholds")
        classif_layout = QtWidgets.QVBoxLayout()

        # Healthy threshold
        self.classif_healthy_label = QtWidgets.QLabel("Healthy : > 5.0 mm")
        self.classif_healthy_label.setStyleSheet("font-size: 13px; font-weight: bold; color: #4caf50;")
        classif_layout.addWidget(self.classif_healthy_label)
        self.classif_healthy_spin = QtWidgets.QDoubleSpinBox()
        self.classif_healthy_spin.setRange(1.0, 20.0)
        self.classif_healthy_spin.setValue(5.0)
        self.classif_healthy_spin.setSingleStep(0.5)
        self.classif_healthy_spin.setSuffix(" mm")
        self.classif_healthy_spin.setKeyboardTracking(False)
        self.classif_healthy_spin.valueChanged.connect(self._on_classif_thresholds_changed)
        classif_layout.addWidget(self.classif_healthy_spin)

        # Border threshold
        self.classif_border_label = QtWidgets.QLabel("Border zone : > 4.0 mm")
        self.classif_border_label.setStyleSheet("font-size: 13px; font-weight: bold; color: #f9a825;")
        classif_layout.addWidget(self.classif_border_label)
        self.classif_border_spin = QtWidgets.QDoubleSpinBox()
        self.classif_border_spin.setRange(0.5, 15.0)
        self.classif_border_spin.setValue(4.0)
        self.classif_border_spin.setSingleStep(0.5)
        self.classif_border_spin.setSuffix(" mm")
        self.classif_border_spin.setKeyboardTracking(False)
        self.classif_border_spin.valueChanged.connect(self._on_classif_thresholds_changed)
        classif_layout.addWidget(self.classif_border_spin)

        # Dense scar threshold
        self.classif_dense_label = QtWidgets.QLabel("Dense scar : ≤ 2.0 mm")
        self.classif_dense_label.setStyleSheet("font-size: 13px; font-weight: bold; color: #c62828;")
        classif_layout.addWidget(self.classif_dense_label)
        self.classif_dense_spin = QtWidgets.QDoubleSpinBox()
        self.classif_dense_spin.setRange(0.1, 10.0)
        self.classif_dense_spin.setValue(2.0)
        self.classif_dense_spin.setSingleStep(0.5)
        self.classif_dense_spin.setSuffix(" mm")
        self.classif_dense_spin.setKeyboardTracking(False)
        self.classif_dense_spin.valueChanged.connect(self._on_classif_thresholds_changed)
        classif_layout.addWidget(self.classif_dense_spin)

        classif_hint = QtWidgets.QLabel("Scar = entre Dense et Border")
        classif_hint.setStyleSheet("font-size: 11px; color: #888;")
        classif_layout.addWidget(classif_hint)

        # Checkbox pour ne visualiser que Scar + BZ
        self.scar_bz_checkbox = QtWidgets.QCheckBox("Scar + BZ uniquement")
        self.scar_bz_checkbox.setToolTip(
            "Masquer les zones saines pour ne garder que\n"
            "Scar + Border Zone + Dense Scar")
        self.scar_bz_checkbox.setChecked(False)
        self.scar_bz_checkbox.stateChanged.connect(self._on_scar_bz_toggled)
        classif_layout.addWidget(self.scar_bz_checkbox)

        self.classif_group.setLayout(classif_layout)
        self.classif_group.setVisible(True)
        layout.addWidget(self.classif_group)

        stats_group = QtWidgets.QGroupBox("Statistics")
        stats_layout = QtWidgets.QVBoxLayout()
        
        self.stats_text = QtWidgets.QTextEdit()
        self.stats_text.setReadOnly(True)
        self.stats_text.setMaximumHeight(450)
        self.stats_text.setText("Load patient data to see statistics")
        
        stats_layout.addWidget(self.stats_text)
        
        dashboard_btn = QtWidgets.QPushButton("Open Dashboard Report")
        dashboard_btn.setStyleSheet("""
            QPushButton {
                background-color: #7b1fa2; color: white;
                padding: 10px 16px; border-radius: 6px;
                font-size: 14px; font-weight: bold;
            }
            QPushButton:hover { background-color: #6a1b9a; }
            QPushButton:pressed { background-color: #4a148c; }
        """)
        dashboard_btn.clicked.connect(self._show_dashboard)
        stats_layout.addWidget(dashboard_btn)
        
        stats_group.setLayout(stats_layout)
        layout.addWidget(stats_group)
        
        # Spacer
        layout.addStretch()

        # Enrouler le panneau dans un QScrollArea pour que les sliders
        # dynamiques (r, p, h_min…) restent toujours accessibles quand
        # la fenêtre est petite.
        scroll = QtWidgets.QScrollArea()
        scroll.setWidgetResizable(True)
        scroll.setHorizontalScrollBarPolicy(QtCore.Qt.ScrollBarAlwaysOff)
        scroll.setVerticalScrollBarPolicy(QtCore.Qt.ScrollBarAsNeeded)
        scroll.setFrameShape(QtWidgets.QFrame.NoFrame)
        scroll.setMinimumWidth(340)
        scroll.setMaximumWidth(430)
        panel.setMinimumWidth(0)
        panel.setMaximumWidth(16777215)
        scroll.setWidget(panel)
        return scroll
    
    # =========================================================================
    # GESTION DYNAMIQUE DE LA GRILLE
    # =========================================================================
    
    def _rebuild_view_combos(self):
        """Reconstruit les combos en fonction du nombre de vues courant"""
        
        # Sauvegarder les sélections existantes
        old_selections = {}
        for idx, combo in self.view_combos.items():
            old_selections[idx] = combo.currentData()
        
        # Supprimer les anciens widgets
        while self.views_inner_layout.count():
            item = self.views_inner_layout.takeAt(0)
            if item.layout():
                while item.layout().count():
                    child = item.layout().takeAt(0)
                    if child.widget():
                        child.widget().deleteLater()
            elif item.widget():
                item.widget().deleteLater()
        
        self.view_combos = {}
        
        for i in range(self._n_views):
            row_layout = QtWidgets.QHBoxLayout()
            
            label = QtWidgets.QLabel(f"View {i+1}:")
            label.setStyleSheet("font-size: 13px; font-weight: 600;")
            row_layout.addWidget(label, stretch=0)
            
            combo = QtWidgets.QComboBox()
            combo.setStyleSheet("font-size: 13px;")
            combo.addItem("-- Empty --", None)
            
            for category, metric_ids in self.catalog.categories.items():
                for metric_id in metric_ids:
                    metric = self.catalog.get_metric(metric_id)
                    combo.addItem(f"{category}: {metric.name}", metric_id)
            
            # Restaurer la sélection si elle existait
            if i in old_selections and old_selections[i] is not None:
                for ci in range(combo.count()):
                    if combo.itemData(ci) == old_selections[i]:
                        combo.setCurrentIndex(ci)
                        break
            
            combo.currentIndexChanged.connect(
                lambda idx, view_idx=i: self._on_view_changed(view_idx, idx)
            )
            
            row_layout.addWidget(combo, stretch=1)
            self.view_combos[i] = combo
            self.views_inner_layout.addLayout(row_layout)
    
    def _create_view_buttons_bar(self):
        """Crée un panneau flottant vertical de raccourcis de vue (bas gauche du plotter)."""
        overlay = QtWidgets.QFrame(self.viz_container)
        overlay.setObjectName("view_overlay")
        overlay.setStyleSheet("""
            QFrame#view_overlay {
                background-color: white;
                border: 1px solid #d0d0d0;
                border-radius: 6px;
            }
        """)
        overlay.setMaximumWidth(320)
        lay = QtWidgets.QVBoxLayout(overlay)
        lay.setContentsMargins(4, 6, 4, 6)
        lay.setSpacing(2)

        BTN_STYLE = """
            QPushButton {
                background: transparent; border: none;
                border-radius: 4px; padding: 3px 6px;
                font-size: 11px; font-weight: 600;
                color: #333; min-width: 40px; min-height: 24px;
            }
            QPushButton:hover  { background-color: #1976D2; color: white; }
            QPushButton:pressed { background-color: #0D47A1; color: white; }
        """
        SEP_STYLE = "QFrame { color: #d0d0d0; background: #d0d0d0; max-height: 1px; margin: 2px 0; }"

        lbl = QtWidgets.QLabel("Vue")
        lbl.setAlignment(QtCore.Qt.AlignCenter)
        lbl.setStyleSheet("font-size: 9px; color: #999; background: transparent; border: none; margin: 0;")
        lbl.setMaximumHeight(14)
        lay.addWidget(lbl)

        sep = QtWidgets.QFrame()
        sep.setFrameShape(QtWidgets.QFrame.HLine)
        sep.setStyleSheet(SEP_STYLE)
        lay.addWidget(sep)

        views = [
            ("AP",  "anterior",  "Antéro-postérieure"),
            ("PA",  "posterior", "Postéro-antérieure"),
            ("LAO", "lao",       "Left Anterior Oblique (40°)"),
            ("RAO", "rao",       "Right Anterior Oblique (40°)"),
            ("LL",  "ll",        "Latérale gauche"),
            ("RL",  "rl",        "Latérale droite"),
            ("INF", "inferior",  "Vue inférieure (apex)"),
            ("SUP", "superior",  "Vue supérieure (base)"),
        ]
        for text, vtype, tip in views:
            btn = QtWidgets.QPushButton(text)
            btn.setToolTip(tip)
            btn.setStyleSheet(BTN_STYLE)
            btn.clicked.connect(lambda _c, vt=vtype: self._set_all_views(vt))
            lay.addWidget(btn)

        sep2 = QtWidgets.QFrame()
        sep2.setFrameShape(QtWidgets.QFrame.HLine)
        sep2.setStyleSheet(SEP_STYLE)
        lay.addWidget(sep2)

        reset_btn = QtWidgets.QPushButton("↺")
        reset_btn.setToolTip("Réinitialiser (vue antérieure)")
        reset_btn.setStyleSheet(BTN_STYLE)
        reset_btn.clicked.connect(self._on_reset_camera)
        lay.addWidget(reset_btn)

        # Positionnement initial + event filter pour le suivi de resize
        overlay.adjustSize()
        self._overlay_mgr = _VizOverlayManager(overlay, self.viz_container)
        self._overlay_mgr._reposition(self.viz_container)
        overlay.show()
        overlay.raise_()
        self._view_overlay = overlay

    def _set_all_views(self, view_type: str):
        """Applique une vue caméra à toutes les sous-vues actives."""
        if not self.viz_manager:
            return
        rows, cols = self._grid_shape
        for r in range(rows):
            for c in range(cols):
                self.viz_manager.set_view((r, c), view_type)
        self.plotter.render()

    def _create_plotter(self, shape: tuple):
        """Crée (ou recrée) le QtInteractor avec la shape donnée"""
        
        # Nettoyer l'ancien widget d'orientation
        if self.orientation_widget is not None:
            self.orientation_widget.cleanup()
            self.orientation_widget = None
        
        # Fermer l'ancien plotter proprement
        if self.plotter is not None:
            try:
                self.plotter.close()
            except:
                pass
            self.plotter.setParent(None)
            self.plotter.deleteLater()
            self.plotter = None
        
        self._grid_shape = shape
        
        self.plotter = QtInteractor(parent=self.viz_container, shape=shape)
        self.plotter.set_background("white")
        # Toujours insérer le plotter en position 0 pour que la barre de vues
        # reste systématiquement en bas, même lors d'un recréation.
        self.viz_layout.insertWidget(0, self.plotter)
        
        # Réinitialiser le viz manager
        self.viz_manager = VisualizationManager(self.plotter)
        self.viz_manager.enable_linked_cursor()
        
        # Fournir le base_mesh au viz manager si on a un data_manager
        if self.data_manager and 'lv_epi_dist' in self.data_manager.meshes:
            self.viz_manager._base_mesh = self.data_manager.meshes['lv_epi_dist']
        
        # Afficher des vues vides
        rows, cols = shape
        for r in range(rows):
            for c in range(cols):
                self.plotter.subplot(r, c)
                self.plotter.add_text("Select metric", font_size=16)
        
        # Synchroniser les vues par défaut
        try:
            self.plotter.link_views()
        except:
            pass
        
        # Activer le picking par clic droit pour choisir le site de pacing
        try:
            self.plotter.enable_point_picking(
                callback=self._on_right_click_pacing,
                use_picker=True,
                show_point=True,
                point_size=14,
                color='yellow',
                show_message=False,
                left_clicking=False,   # clic droit uniquement
            )
        except Exception:
            pass
        
        # Créer le widget d'orientation avec buste humain (bas droite, unique)
        try:
            self.orientation_widget = HumanBustOrientationWidget()
            # Récupérer l'interactor VTK du plotter pyvistaqt
            iren = self.plotter.iren
            if hasattr(iren, 'interactor'):
                vtk_iren = iren.interactor
            else:
                vtk_iren = iren
            self.orientation_widget.setup(vtk_iren)
        except Exception:
            self.orientation_widget = None
        
        # Mettre à jour les labels
        if hasattr(self, 'grid_label'):
            self.grid_label.setText(f"({rows}×{cols})")
        if hasattr(self, 'view_count_label'):
            self.view_count_label.setText(str(self._n_views))

        # Réinitialiser le tracking caméra (nouvelle grille = nouvelles vues vierges)
        self._initialized_views = set()

        # Remonter l'overlay de vues au premier plan
        if hasattr(self, '_view_overlay') and self._view_overlay:
            self._view_overlay.raise_()
    
    def _on_add_view(self):
        """Ajoute une vue (bouton +)"""
        if self._n_views < self.MAX_VIEWS:
            self._set_view_count(self._n_views + 1)
    
    def _on_remove_view(self):
        """Retire une vue (bouton −)"""
        if self._n_views > 1:
            self._set_view_count(self._n_views - 1)
    
    def _set_view_count(self, n: int):
        """Change le nombre de vues et reconstruit la grille"""
        self._n_views = n
        new_shape = _optimal_grid(n)
        
        # Recréer combos et plotter
        self._rebuild_view_combos()
        self._create_plotter(new_shape)
        
        # Re-appliquer les métriques qui étaient sélectionnées
        self._reapply_all_views()
    
    def _auto_adjust_views(self):
        """Ajuste automatiquement le nombre de vues au nombre de métriques sélectionnées"""
        active_count = sum(
            1 for combo in self.view_combos.values()
            if combo.currentData() is not None
        )
        
        if active_count == 0:
            active_count = 1
        
        self._set_view_count(active_count)
    
    def _reapply_all_views(self):
        """Réapplique toutes les métriques sélectionnées dans les combos"""
        if not self.data_manager:
            return
        
        for view_idx, combo in self.view_combos.items():
            combo_idx = combo.currentIndex()
            if combo_idx > 0:
                self._on_view_changed(view_idx, combo_idx)
        
        self.plotter.render()
    
    def _create_menu_bar(self):
        """Crée la barre de menu"""
        
        menubar = self.menuBar()
        
        # Menu File
        file_menu = menubar.addMenu("File")
        
        load_action = QtGui.QAction("Load Patient...", self)
        load_action.triggered.connect(self._on_load_patient)
        file_menu.addAction(load_action)
        
        file_menu.addSeparator()
        
        quit_action = QtGui.QAction("Quit", self)
        quit_action.triggered.connect(self.close)
        file_menu.addAction(quit_action)
        
        # Menu View
        view_menu = menubar.addMenu("View")
        
        add_view_action = QtGui.QAction("Add View (+)", self)
        add_view_action.triggered.connect(self._on_add_view)
        view_menu.addAction(add_view_action)
        
        remove_view_action = QtGui.QAction("Remove View (−)", self)
        remove_view_action.triggered.connect(self._on_remove_view)
        view_menu.addAction(remove_view_action)        
        # Menu Analysis
        analysis_menu = menubar.addMenu("Analysis")
        
        dashboard_action = QtGui.QAction("Open Dashboard...", self)
        dashboard_action.triggered.connect(self._show_dashboard)
        analysis_menu.addAction(dashboard_action)    
    # === CALLBACKS ===
    
    def _on_load_patient(self):
        """Callback : Charger un patient"""
        
        folder = QtWidgets.QFileDialog.getExistingDirectory(
            self,
            "Select Patient Data Folder",
            "",
            QtWidgets.QFileDialog.ShowDirsOnly
        )
        
        if folder:
            self.load_patient(folder)
    
    def load_patient(self, data_path: str):
        """Charge les données d'un patient"""
        
        try:
            self.data_path = data_path
            self.data_manager = DataManager(data_path)
            
            # Charger les données
            success = self.data_manager.load_patient_data()
            
            if not success:
                QtWidgets.QMessageBox.warning(
                    self,
                    "Load Error",
                    "No VTK files found in the selected folder"
                )
                return
            
            # Fournir le mesh EPI de base au viz manager (pour ghost overlays)
            if self.viz_manager and 'lv_epi_dist' in self.data_manager.meshes:
                self.viz_manager._base_mesh = self.data_manager.meshes['lv_epi_dist']
            
            # Mettre à jour l'interface
            self.patient_label.setText(f"Loaded: {data_path.split('/')[-1]}")
            
            # Calculer et afficher les statistiques
            self._update_statistics()
            
            # Message de succès
            QtWidgets.QMessageBox.information(
                self,
                "Success",
                f"Patient data loaded successfully!\n{len(self.data_manager.meshes)} meshes found."
            )
            
        except Exception as e:
            QtWidgets.QMessageBox.critical(
                self,
                "Error",
                f"Failed to load patient data:\n{str(e)}"
            )
    
    def _update_statistics(self):
        """Met à jour l'affichage des statistiques"""
        
        if not self.data_manager:
            return
        
        stats = self.data_manager.get_statistics()
        
        text = "\u2550" * 32 + "\n"
        text += "   VENTRICULAR REPORT\n"
        text += "\u2550" * 32 + "\n\n"
        
        # --- Geometry ---
        text += "\u25b8 GEOMETRY\n"
        if 'LV_Surface_cm2' in stats:
            text += f"  Surface:    {stats['LV_Surface_cm2']:>8.1f} cm\u00b2\n"
        if 'LV_Volume_mL' in stats:
            text += f"  Volume:     {stats['LV_Volume_mL']:>8.1f} mL\n"
        if 'LV_Long_Axis_mm' in stats:
            text += f"  Long Axis:  {stats['LV_Long_Axis_mm']:>8.1f} mm\n"
        if 'LV_Sphericity' in stats:
            text += f"  Sphericity: {stats['LV_Sphericity']:>8.3f}\n"
        if 'LV_N_Points' in stats:
            text += f"  Mesh pts:   {stats['LV_N_Points']:>8,}\n"
        text += "\n"
        
        # --- Thickness ---
        text += "\u25b8 WALL THICKNESS\n"
        if 'T_mean' in stats:
            text += f"  Mean \u00b1 SD:  {stats['T_mean']:.2f} \u00b1 {stats['T_std']:.2f} mm\n"
            text += f"  Median:     {stats['T_median']:>8.2f} mm\n"
            text += f"  Range:     [{stats['T_min']:.1f} \u2013 {stats['T_max']:.1f}] mm\n"
            text += f"  IQR:        {stats.get('T_IQR', 0):>8.2f} mm\n"
            text += f"  Scar (<1):  {stats['Scar_Pct']:>7.1f}%\n"
            text += f"  Border(1-5):{stats['Border_Zone_Pct']:>7.1f}%\n"
            text += f"  Healthy(>5):{stats['Healthy_Pct']:>7.1f}%\n"
        text += "\n"
        
        # --- Arrhythmia ---
        text += "\u25b8 ARRHYTHMIA SUBSTRATE\n"
        if 'DZ_Area_Pct' in stats:
            text += f"  DZ extent:  {stats['DZ_Area_Pct']:>7.1f}%\n"
            text += f"  Max \u03c1:      {stats.get('Rho_Max', 0):>8.4f}\n"
        if 'Isthmus_N_Points' in stats:
            text += f"  Isthmus pts:{stats['Isthmus_N_Points']:>8d}\n"
        if 'Channel_Pct' in stats:
            text += f"  Channels:   {stats['Channel_Pct']:>7.1f}%\n"
        text += "\n"
        
        # --- Scar ---
        text += "\u25b8 SCAR BURDEN (Utah)\n"
        if 'Scar_Burden_Pct' in stats:
            text += f"  Scar/Total: {stats['Scar_Burden_Pct']:>7.1f}%\n"
            text += f"  Dense/Total:{stats.get('Dense_Burden_Pct', 0):>7.1f}%\n"
            text += f"  Utah Grade: {stats.get('Utah_Grade', 'N/A')}\n"
            text += f"  Scar Area:  {stats.get('Scar_Area_CT_cm2', 0):>7.2f} cm\u00b2\n"
            text += f"  Dense Area: {stats.get('Dense_Area_CT_cm2', 0):>7.2f} cm\u00b2\n"
        elif 'Dense_Scar_Surface_cm2' in stats:
            text += f"  Surface:    {stats['Dense_Scar_Surface_cm2']:>8.2f} cm\u00b2\n"
            if 'Dense_Scar_Burden_Pct' in stats:
                text += f"  Burden:     {stats['Dense_Scar_Burden_Pct']:>7.1f}%\n"
            if 'Dense_Scar_Tissue_Volume_cm3' in stats:
                text += f"  Tissue Vol: {stats['Dense_Scar_Tissue_Volume_cm3']:>6.2f} cm\u00b3\n"
        text += "\n"
        
        # --- Channelness ---
        text += "\u25b8 CHANNELNESS\n"
        if 'Channelness_Mean' in stats:
            text += f"  Mean:       {stats['Channelness_Mean']:>8.4f}\n"
            text += f"  Max:        {stats['Channelness_Max']:>8.4f}\n"
            text += f"  High (>0.5):{stats.get('Channelness_High_Pct', 0):>7.1f}%\n"
        else:
            text += "  (not computed)\n"
        text += "\n"
        
        # --- Border Zone Cedilnik ---
        text += "\u25b8 BORDER ZONE (Cedilnik)\n"
        if 'Cedilnik_BZ_Pct' in stats:
            text += f"  BZ extent:  {stats['Cedilnik_BZ_Pct']:>7.1f}%\n"
            text += f"  BZ area:    {stats.get('Cedilnik_BZ_Area_cm2', 0):>7.2f} cm\u00b2\n"
            text += f"  Radius:     {stats.get('Cedilnik_Radius_mm', 2.0):>7.1f} mm\n"
        else:
            text += "  (not computed)\n"
        text += "\n"
        
        # --- Transmurality ---
        text += "\u25b8 TRANSMURALITY\n"
        if 'Trans_Mean' in stats:
            text += f"  Mean:       {stats['Trans_Mean']:>8.1f}%\n"
            if 'Trans_Transmural_Pct' in stats:
                text += f"  Transmural: {stats['Trans_Transmural_Pct']:>7.1f}%\n"
                text += f"  Midmural:   {stats['Trans_Midmural_Pct']:>7.1f}%\n"
        text += "\n"
        
        text += "\u2550" * 32 + "\n"
        text += "  Open Dashboard for full\n"
        text += "  report & charts\n"
        text += "\u2550" * 32
        
        self.stats_text.setText(text)
    
    def _show_dashboard(self):
        """Ouvre la fen\u00eatre de dashboard avec graphiques complets"""
        if not self.data_manager:
            QtWidgets.QMessageBox.warning(
                self, "No Data",
                "Please load patient data first")
            return
        
        if not _HAS_DASHBOARD:
            QtWidgets.QMessageBox.warning(
                self, "Missing Dependency",
                "matplotlib is required for the dashboard.\n"
                "Install with: pip install matplotlib")
            return
        
        try:
            dashboard = DashboardWindow(self.data_manager, parent=self)
            dashboard.exec_()
        except Exception as e:
            QtWidgets.QMessageBox.critical(
                self, "Dashboard Error",
                f"Failed to generate dashboard:\n{str(e)}")
    
    def _on_view_changed(self, view_idx: int, combo_idx: int):
        """Callback : Une vue a changé de métrique — ne touche PAS aux autres vues"""
        
        if not self.data_manager:
            QtWidgets.QMessageBox.warning(
                self,
                "No Data",
                "Please load patient data first"
            )
            # Reset combo sans re-déclencher le signal
            combo = self.view_combos[view_idx]
            combo.blockSignals(True)
            combo.setCurrentIndex(0)
            combo.blockSignals(False)
            return
        
        # Récupérer la métrique sélectionnée
        combo = self.view_combos[view_idx]
        metric_id = combo.itemData(combo_idx)
        
        # Position dans la grille adaptive
        rows, cols = self._grid_shape
        i = view_idx // cols
        j = view_idx % cols
        
        # Vérifier que l'index est dans la grille
        if i >= rows or j >= cols:
            return
        
        if metric_id is None:
            # Vue vide
            self.viz_manager.render_empty((i, j), "Select metric")
        else:
            # Récupérer la définition de la métrique
            metric_def = self.catalog.get_metric(metric_id)
            
            # Calculer/récupérer le mesh
            compute_method = getattr(self.data_manager, metric_def.compute_func)
            
            # Gestion des paramètres pour les simulations
            if metric_id == 'simulation':
                decay = self.scar_decay_slider.value() / 2.0
                custom_pt = getattr(self, '_custom_pacing_point', None)
                if custom_pt is not None:
                    mesh = compute_method(pacing_site='CUSTOM', scar_decay=decay, custom_point=custom_pt)
                else:
                    mesh = compute_method(pacing_site='SR', scar_decay=decay)
            elif metric_id == 'channels':
                border_radius = self.border_slider.value()
                mesh = compute_method(border_radius=border_radius)
            elif metric_id == 'border_zone_cedilnik':
                dilation_radius = self.cedilnik_slider.value() / 2.0
                mesh = compute_method(dilation_radius=dilation_radius)
            elif metric_id == 'channelness':
                p_val = self.channelness_slider.value() / 10.0
                mesh = compute_method(p=p_val)
            elif metric_id == 'anatomical_channelness':
                h_min     = self.anat_ch_hmin_slider.value() / 10.0
                max_width = self.anat_ch_maxw_slider.value() / 10.0
                if 'anatomical_channelness' in self.data_manager.computed_metrics:
                    del self.data_manager.computed_metrics['anatomical_channelness']
                mesh = compute_method(h_min=h_min, max_width=max_width)
            elif metric_id == 'scar_burden':
                thresh_h = self.classif_healthy_spin.value()
                thresh_b = self.classif_border_spin.value()
                thresh_d = self.classif_dense_spin.value()
                if 'scar_burden' in self.data_manager.computed_metrics:
                    del self.data_manager.computed_metrics['scar_burden']
                mesh = compute_method(thresh_healthy=thresh_h, thresh_border=thresh_b, thresh_dense=thresh_d)
            elif metric_id == 'combined_score':
                # Ouvrir le dialogue de configuration
                prev_config = self._combined_configs.get(view_idx)
                dlg = CombinedScoreDialog(parent=self, previous_config=prev_config)
                if dlg.exec_() != QtWidgets.QDialog.Accepted:
                    # Annulé — remettre combo à l'ancienne valeur
                    combo_widget = self.view_combos[view_idx]
                    combo_widget.blockSignals(True)
                    combo_widget.setCurrentIndex(0)
                    combo_widget.blockSignals(False)
                    return
                cfg = dlg.get_config()
                self._combined_configs[view_idx] = cfg
                mesh = compute_method(config=cfg)
            else:
                mesh = compute_method()
            
            if mesh is None:
                self.viz_manager.render_empty((i, j), f"Data unavailable:\n{metric_def.name}")
            else:
                # Appliquer le filtre Scar+BZ si activé
                mesh = self._apply_scar_bz_filter(mesh)

                # Conserver la caméra si le subplot avait déjà un métrique
                is_fresh = (i, j) not in self._initialized_views
                saved_cam = None
                if not is_fresh:
                    try:
                        self.plotter.subplot(i, j)
                        saved_cam = list(self.plotter.camera_position)
                    except Exception:
                        pass

                # Paramètres spécifiques par métrique
                extra_kwargs = {}
                if metric_id == 'deceleration':
                    extra_kwargs['dz_threshold'] = self.dz_slider.value() / 100.0
                elif metric_id == 'laplacian':
                    extra_kwargs['lap_threshold'] = self.lap_thresh_spin.value()
                    extra_kwargs['lap_min_area'] = self.lap_area_spin.value()
                
                # Afficher la métrique — ne touche qu'à ce subplot
                self.viz_manager.render_metric((i, j), mesh, metric_def, **extra_kwargs)

                # Caméra : vue antérieure si premier chargement, sinon restaurer
                if is_fresh:
                    self._initialized_views.add((i, j))
                    self.viz_manager.set_view((i, j), 'anterior')
                elif saved_cam is not None:
                    try:
                        self.plotter.subplot(i, j)
                        self.plotter.camera_position = saved_cam
                    except Exception:
                        pass

        # Afficher/masquer le slider DZ selon si une vue affiche la métrique deceleration
        has_dz = any(
            combo.currentData() == 'deceleration'
            for combo in self.view_combos.values()
        )
        self.dz_group.setVisible(has_dz)

        # Afficher/masquer le slider Épaisseur
        has_thickness = any(
            combo.currentData() in ('thickness', 'parietal_thickness')
            for combo in self.view_combos.values()
        )
        self.thickness_group.setVisible(has_thickness)

        # Afficher/masquer le slider Border Zone selon si une vue affiche channels
        has_channels = any(
            combo.currentData() == 'channels'
            for combo in self.view_combos.values()
        )
        self.border_group.setVisible(has_channels)
        
        # Afficher/masquer le slider Laplacian
        has_lap = any(
            combo.currentData() == 'laplacian'
            for combo in self.view_combos.values()
        )
        self.lap_group.setVisible(has_lap)
        
        # Afficher/masquer le slider Cedilnik BZ
        has_cedilnik = any(
            combo.currentData() == 'border_zone_cedilnik'
            for combo in self.view_combos.values()
        )
        self.cedilnik_group.setVisible(has_cedilnik)
        
        # Afficher/masquer le slider Channelness
        has_channelness = any(
            combo.currentData() == 'channelness'
            for combo in self.view_combos.values()
        )
        self.channelness_group.setVisible(has_channelness)

        # Afficher/masquer le slider Anatomical Channelness
        has_anat_ch = any(
            combo.currentData() == 'anatomical_channelness'
            for combo in self.view_combos.values()
        )
        self.anat_ch_group.setVisible(has_anat_ch)

        # Afficher/masquer le groupe Simulation Activation
        has_activation = any(
            combo.currentData() == 'simulation'
            for combo in self.view_combos.values()
        )
        self.activation_group.setVisible(has_activation)

        # Le groupe classification est toujours visible (utilisé pour le filtre scar+BZ)

        # Re-rendre (tout le plotter, mais les autres subplots gardent leur contenu)
        self.plotter.render()
    
    def _on_thickness_slider_changed(self, value):
        """Callback : Le slider épaisseur a changé — met à jour clim en temps réel (sans clignotement)."""
        clim_max = value / 10.0
        step = clim_max / 6.0
        self.thickness_label.setText(f"Max : {clim_max:.1f} mm  (6 zones × {step:.2f} mm)")

        if not self.data_manager or not self.viz_manager:
            return

        try:
            rows, cols = self._grid_shape
            for view_idx, combo in self.view_combos.items():
                metric_id = combo.currentData()
                if metric_id in ('thickness', 'parietal_thickness'):
                    i = view_idx // cols
                    j = view_idx % cols
                    if i >= rows or j >= cols:
                        continue
                    # Mettre à jour le clim sans re-rendre (juste la barre de couleur change)
                    self.viz_manager.set_clim_on_actors((i, j), (0.0, clim_max))
            # Un seul rendu en fin pour mettre à jour l'affichage
            self.plotter.render()
        except Exception:
            pass

    def _on_dz_threshold_changed(self, value):
        """Callback : Le slider DZ a changé — re-rendre toutes les vues deceleration"""
        threshold = value / 100.0
        self.dz_label.setText(f"Threshold: {threshold:.2f}")
        
        if not self.data_manager:
            return
        
        try:
            rows, cols = self._grid_shape
            for view_idx, combo in self.view_combos.items():
                metric_id = combo.currentData()
                if metric_id == 'deceleration':
                    i = view_idx // cols
                    j = view_idx % cols
                    
                    if i >= rows or j >= cols:
                        continue
                    
                    metric_def = self.catalog.get_metric(metric_id)
                    mesh = self.data_manager.compute_deceleration_zones(dz_threshold=threshold)
                    
                    if mesh is None:
                        continue
                    
                    self.viz_manager.render_metric((i, j), mesh, metric_def, dz_threshold=threshold)
            
            self.plotter.render()
        
        except Exception:
            pass

    def _on_border_radius_changed(self, value):
        """Callback : Le slider border radius a changé — re-rendre toutes les vues channels"""
        self.border_label.setText(f"Radius: {value:.0f} mm")
        
        if not self.data_manager:
            return
        
        try:
            rows, cols = self._grid_shape
            for view_idx, combo in self.view_combos.items():
                metric_id = combo.currentData()
                if metric_id == 'channels':
                    i = view_idx // cols
                    j = view_idx % cols
                    
                    if i >= rows or j >= cols:
                        continue
                    
                    metric_def = self.catalog.get_metric(metric_id)
                    mesh = self.data_manager.compute_channels(border_radius=float(value))
                    
                    if mesh is None:
                        continue
                    
                    self.viz_manager.render_metric((i, j), mesh, metric_def)
            
            self.plotter.render()
        
        except Exception:
            pass
    
    def _on_lap_threshold_changed(self, value=None):
        """Callback : une spinbox Laplacian a changé — re-rendre toutes les vues laplacian"""
        threshold = self.lap_thresh_spin.value()
        min_area  = self.lap_area_spin.value()

        if not self.data_manager:
            return

        try:
            rows, cols = self._grid_shape
            for view_idx, combo in self.view_combos.items():
                metric_id = combo.currentData()
                if metric_id != 'laplacian':
                    continue
                i = view_idx // cols
                j = view_idx % cols
                if i >= rows or j >= cols:
                    continue
                metric_def = self.catalog.get_metric(metric_id)
                mesh = self.data_manager.compute_laplacian()
                if mesh is None:
                    continue
                mesh = self._apply_scar_bz_filter(mesh)
                # Conserver la caméra
                saved_cam = None
                try:
                    self.plotter.subplot(i, j)
                    saved_cam = list(self.plotter.camera_position)
                except Exception:
                    pass
                self.viz_manager.render_metric((i, j), mesh, metric_def,
                                               lap_threshold=threshold,
                                               lap_min_area=min_area)
                if saved_cam is not None:
                    try:
                        self.plotter.subplot(i, j)
                        self.plotter.camera_position = saved_cam
                    except Exception:
                        pass
            self.plotter.render()
        except Exception:
            pass
    
    def _on_cedilnik_radius_changed(self, value):
        """Callback : Le slider Cedilnik a changé — re-rendre toutes les vues border_zone_cedilnik"""
        radius = value / 2.0
        self.cedilnik_label.setText(f"Radius: {radius:.1f} mm")
        
        if not self.data_manager:
            return
        
        try:
            rows, cols = self._grid_shape
            for view_idx, combo in self.view_combos.items():
                metric_id = combo.currentData()
                if metric_id == 'border_zone_cedilnik':
                    i = view_idx // cols
                    j = view_idx % cols
                    
                    if i >= rows or j >= cols:
                        continue
                    
                    metric_def = self.catalog.get_metric(metric_id)
                    # Forcer le recalcul avec le nouveau rayon
                    if 'border_zone_cedilnik' in self.data_manager.computed_metrics:
                        del self.data_manager.computed_metrics['border_zone_cedilnik']
                    mesh = self.data_manager.compute_border_zone_cedilnik(dilation_radius=radius)
                    
                    if mesh is None:
                        continue
                    
                    self.viz_manager.render_metric((i, j), mesh, metric_def)
            
            self.plotter.render()
        
        except Exception:
            pass
    
    def _on_channelness_params_changed(self, _value=None):
        """Callback : p (inflexion mm) ou r (raideur mm⁻¹) a changé — recalcul."""
        p_val = self.channelness_slider.value() / 10.0
        r_val = self.channelness_r_slider.value() / 10.0
        self.channelness_label.setText(f"Seuil p : {p_val:.1f} mm")
        self.channelness_r_label.setText(f"Raideur r : {r_val:.1f} mm\u207b\u00b9")

        if not self.data_manager:
            return

        try:
            rows, cols = self._grid_shape
            for view_idx, combo in self.view_combos.items():
                metric_id = combo.currentData()
                if metric_id == 'channelness':
                    i = view_idx // cols
                    j = view_idx % cols

                    if i >= rows or j >= cols:
                        continue

                    metric_def = self.catalog.get_metric(metric_id)
                    if 'channelness' in self.data_manager.computed_metrics:
                        del self.data_manager.computed_metrics['channelness']
                    mesh = self.data_manager.compute_channelness(p=p_val, r=r_val)

                    if mesh is None:
                        continue

                    self.viz_manager.render_metric((i, j), mesh, metric_def)

            self.plotter.render()

        except Exception:
            pass

    def _on_scar_decay_changed(self, value):
        """Callback : le slider decay cicatrice a changé — invalide substrat + simulation."""
        decay = value / 2.0
        self.scar_decay_label.setText(f"Decay cicatrice : {decay:.1f} mm")

        if not self.data_manager:
            return

        try:
            # Invalider le substrat et tous les caches simulation
            self.data_manager.computed_metrics.pop('_substrate', None)
            self.data_manager.computed_metrics.pop('simulation', None)

            rows, cols = self._grid_shape
            for view_idx, combo in self.view_combos.items():
                metric_id = combo.currentData() or ''
                if metric_id != 'simulation':
                    continue
                i = view_idx // cols
                j = view_idx % cols
                if i >= rows or j >= cols:
                    continue
                metric_def = self.catalog.get_metric(metric_id)
                compute_method = getattr(self.data_manager, metric_def.compute_func)
                custom_pt = getattr(self, '_custom_pacing_point', None)
                if custom_pt is not None:
                    mesh = compute_method(pacing_site='CUSTOM', scar_decay=decay, custom_point=custom_pt)
                else:
                    mesh = compute_method(pacing_site='SR', scar_decay=decay)
                if mesh is None:
                    continue
                self.viz_manager.render_metric((i, j), mesh, metric_def)

            self.plotter.render()
        except Exception:
            pass

    def _on_at_scale_changed(self, value):
        """Callback : le slider activation time scale a changé — met à jour clim."""
        at_scale = float(value)
        self.at_scale_label.setText(f"Échelle latence : {at_scale:.0f} ms")

        if not self.data_manager:
            return

        try:
            # Mettre à jour la clim sur les vues affichant la simulation
            rows, cols = self._grid_shape
            for view_idx, combo in self.view_combos.items():
                metric_id = combo.currentData() or ''
                if metric_id != 'simulation':
                    continue
                i = view_idx // cols
                j = view_idx % cols
                if i >= rows or j >= cols:
                    continue
                # Mettre à jour clim = [0, at_scale]
                self.viz_manager.set_clim_on_actors((i, j), clim=(0.0, at_scale))

            self.plotter.render()
        except Exception:
            pass

    def _on_anat_ch_changed(self):
        """Callback: h_min ou max_width a changé — re-calcul anatomical channelness."""
        h_min     = self.anat_ch_hmin_slider.value() / 10.0
        max_width = self.anat_ch_maxw_slider.value() / 10.0
        self.anat_ch_hmin_label.setText(f"h_min: {h_min:.1f} mm")
        self.anat_ch_maxw_label.setText(f"max W: {max_width:.1f} mm")

        if not self.data_manager:
            return

        try:
            rows, cols = self._grid_shape
            for view_idx, combo in self.view_combos.items():
                if combo.currentData() != 'anatomical_channelness':
                    continue
                i = view_idx // cols
                j = view_idx % cols
                if i >= rows or j >= cols:
                    continue
                metric_def = self.catalog.get_metric('anatomical_channelness')
                if 'anatomical_channelness' in self.data_manager.computed_metrics:
                    del self.data_manager.computed_metrics['anatomical_channelness']
                mesh = self.data_manager.compute_anatomical_channelness(
                    h_min=h_min, max_width=max_width)
                if mesh is None:
                    continue
                self.viz_manager.render_metric((i, j), mesh, metric_def)
            self.plotter.render()
        except Exception:
            pass

    # ------------------------------------------------------------------
    # Filtre Scar + BZ
    # ------------------------------------------------------------------
    def _apply_scar_bz_filter(self, mesh):
        """Si la case 'Scar + BZ uniquement' est cochée, ne garde que les
        cellules dont l'épaisseur de paroi ≤ seuil healthy.
        Retourne le mesh filtré (ou le mesh original si désactivé)."""
        if not self.scar_bz_checkbox.isChecked():
            return mesh

        thresh_h = self.classif_healthy_spin.value()

        # Il faut le Wall_Thickness sur ce mesh.  Deux cas :
        #  - le mesh l'a déjà (scar_burden, thickness, …)
        #  - sinon on le transfère depuis wall_thickness via correspondance géo
        try:
            if 'Wall_Thickness' not in mesh.point_data:
                wt_mesh = self.data_manager.computed_metrics.get('wall_thickness')
                if wt_mesh is None:
                    self.data_manager.compute_wall_thickness()
                    wt_mesh = self.data_manager.computed_metrics.get('wall_thickness')
                if wt_mesh is None:
                    return mesh
                from scipy.spatial import cKDTree
                tree = cKDTree(wt_mesh.points)
                _, idx = tree.query(mesh.points, k=1)
                mesh['_WT_filter'] = wt_mesh['Wall_Thickness'][idx]
                wt_key = '_WT_filter'
            else:
                wt_key = 'Wall_Thickness'

            # Seuil : garder indices de points ≤ thresh_healthy
            mask = mesh[wt_key] <= thresh_h
            if not np.any(mask):
                return mesh  # tout est sain → on renvoie quand même le mesh complet
            point_ids = np.where(mask)[0]
            filtered = mesh.extract_points(point_ids, adjacent_cells=True)
            # Nettoyer le champ temporaire
            if '_WT_filter' in filtered.point_data:
                del filtered.point_data['_WT_filter']
            if filtered.n_points == 0:
                return mesh
            return filtered
        except Exception:
            return mesh

    def _on_scar_bz_toggled(self, _state=None):
        """Callback : la case Scar+BZ uniquement a basculé → re-rendre toutes les vues."""
        self._refresh_all_views()

    def _refresh_all_views(self):
        """Re-render toutes les vues actives (utile après changement de filtre global)."""
        if not self.data_manager or not self.viz_manager:
            return
        for view_idx, combo in self.view_combos.items():
            ci = combo.currentIndex()
            if combo.currentData() is not None:
                # Simuler un changement de vue pour forcer le re-rendu complet
                self._on_view_changed(view_idx, ci)

    def _on_classif_thresholds_changed(self, _value=None):
        """Callback : un des seuils de classification a changé — re-rendre scar_burden
        et, si le filtre Scar+BZ est actif, toutes les vues."""
        if not self.data_manager:
            return

        thresh_h = self.classif_healthy_spin.value()
        thresh_b = self.classif_border_spin.value()
        thresh_d = self.classif_dense_spin.value()

        # Mettre à jour les labels
        self.classif_healthy_label.setText(f"Healthy : > {thresh_h:.1f} mm")
        self.classif_border_label.setText(f"Border zone : > {thresh_b:.1f} mm")
        self.classif_dense_label.setText(f"Dense scar : \u2264 {thresh_d:.1f} mm")

        # Si le filtre scar+BZ est actif, toutes les vues doivent être rafraîchies
        if self.scar_bz_checkbox.isChecked():
            self.data_manager.computed_metrics.pop('scar_burden', None)
            self._refresh_all_views()
            return

        try:
            # Invalider le cache
            self.data_manager.computed_metrics.pop('scar_burden', None)

            rows, cols = self._grid_shape
            for view_idx, combo in self.view_combos.items():
                metric_id = combo.currentData()
                if metric_id != 'scar_burden':
                    continue
                i = view_idx // cols
                j = view_idx % cols
                if i >= rows or j >= cols:
                    continue
                metric_def = self.catalog.get_metric(metric_id)
                mesh = self.data_manager.compute_scar_burden(
                    thresh_healthy=thresh_h,
                    thresh_border=thresh_b,
                    thresh_dense=thresh_d)
                if mesh is None:
                    continue
                self.viz_manager.render_metric((i, j), mesh, metric_def)
            self.plotter.render()
        except Exception:
            pass

    def _on_right_click_pacing(self, point, picker=None):
        """Callback : clic droit sur le mesh pour choisir un site de pacing personnalisé."""
        import numpy as np
        if point is None or not self.data_manager:
            return
        self._custom_pacing_point = np.array(point)
        self.pacing_info_label.setText(
            f"Site de pacing : ({point[0]:.1f}, {point[1]:.1f}, {point[2]:.1f})")

        # Invalider et relancer la simulation
        self.data_manager.computed_metrics.pop('simulation', None)
        self.data_manager.computed_metrics.pop('_substrate', None)
        decay = self.scar_decay_slider.value() / 2.0

        try:
            rows, cols = self._grid_shape
            for view_idx, combo in self.view_combos.items():
                metric_id = combo.currentData() or ''
                if metric_id != 'simulation':
                    continue
                i = view_idx // cols
                j = view_idx % cols
                if i >= rows or j >= cols:
                    continue
                metric_def = self.catalog.get_metric(metric_id)
                compute_method = getattr(self.data_manager, metric_def.compute_func)
                mesh = compute_method(pacing_site='CUSTOM', scar_decay=decay,
                                      custom_point=self._custom_pacing_point)
                if mesh is None:
                    continue
                self.viz_manager.render_metric((i, j), mesh, metric_def)
            self.plotter.render()
        except Exception:
            pass

    def _toggle_left_panel(self):
        """Masque ou affiche le panneau gauche avec une animation fluide."""
        # Stopper toute animation en cours
        if self._panel_anim is not None:
            self._panel_anim.stop()

        if self._panel_collapsed:
            # Déplier
            self.left_panel.setMinimumWidth(0)
            self.left_panel.setMaximumWidth(16777215)   # reset contrainte max
            anim = QPropertyAnimation(self.left_panel, b"maximumWidth")
            anim.setDuration(260)
            anim.setEasingCurve(QEasingCurve.InOutCubic)
            anim.setStartValue(0)
            anim.setEndValue(420)
            anim.finished.connect(lambda: (
                self.left_panel.setMinimumWidth(340),
                self.left_panel.setMaximumWidth(420),
            ))
            self._panel_anim = anim
            anim.start()
            self._panel_toggle_btn.setText("◄")
            self._panel_collapsed = False
        else:
            # Replier
            self.left_panel.setMinimumWidth(0)
            anim = QPropertyAnimation(self.left_panel, b"maximumWidth")
            anim.setDuration(260)
            anim.setEasingCurve(QEasingCurve.InOutCubic)
            anim.setStartValue(self.left_panel.width())
            anim.setEndValue(0)
            self._panel_anim = anim
            anim.start()
            self._panel_toggle_btn.setText("►")
            self._panel_collapsed = True

    def _on_reset_camera(self):
        """Callback : Reset caméra sur toutes les vues (vue antérieure)"""
        
        rows, cols = self._grid_shape
        for r in range(rows):
            for c in range(cols):
                self.viz_manager.set_view((r, c), 'anterior')
        
        self.plotter.render()


def launch_app(data_path: str = None):
    """Lance l'application"""
    
    app = QtWidgets.QApplication.instance()
    if not app:
        app = QtWidgets.QApplication(sys.argv)
    
    window = LVExplorerApp(data_path)
    window.show()
    
    sys.exit(app.exec_())


if __name__ == "__main__":
    launch_app()