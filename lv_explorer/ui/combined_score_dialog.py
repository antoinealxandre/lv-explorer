"""
Dialogue de configuration du score combiné.
Permet de sélectionner les métriques et leurs poids.
"""

from PyQt5 import QtWidgets, QtCore


class CombinedScoreDialog(QtWidgets.QDialog):
    """
    Dialogue pour configurer le score composite :
      - Cases à cocher par métrique
      - Slider de pondération (0 – 100 %)
    """

    # Ordre d'affichage + méta-infos pour l'UI
    METRIC_ROWS = [
        ('narrow_channel',   'Canal étroit (WT 1-4mm)',
         'Zone WT entre 1 et 4mm — couloir de conduction critique\nCedilnik 2018, de Bakker 1993'),
        ('channelness',      'Channelness (Cedilnik)',
         'Score de channelness par simulation Eikonal\nCedilnik, EP-Europace 2018'),
        ('ciaccio',          'Gradient (Ciaccio)',
         'Changement brutal d\'épaisseur (ρ = δT / r·T)\nCiaccio 2024'),
        ('isthmus',          'Isthme prédit (Takigawa)',
         'Probabilité d\'isthme de TV — score multi-critères\nTakigawa 2019 JACC'),
        ('border_zone',      'Border zone (Cedilnik)',
         'Zone bordante dilatée (intersection scar+healthy)\nCedilnik 2018'),
        ('deceleration',     'Zones de décélération',
         'Zones ρ > seuil — ralentissement du front d\'onde\nCiaccio 2024'),
        ('local_entropy',    'Entropie locale',
         'Irrégularité morphologique locale — hétérogénéité\nHsia 2006'),
        ('tri',              'TRI (Rugosité terrain)',
         'Terrain Ruggedness Index — variation d\'épaisseur voisinage\nRiley 1999'),
        ('cv_slow',          'Conduction lente (1 - CV)',
         'Inverse de la vitesse de conduction — zones à risque de bloc\nCedilnik 2018'),
        ('scar_proximity',   'Proximité scar dense',
         'Distance normalisée au scar dense (T ≤ 2mm)\nStevenson 1989 Circ'),
        ('activation_rv_late', 'Activation tardive — Pacing VD',
         'Zones à activation retardée lors d\'un pacing VD (normalisé 0→1)\nSimulation Eikonal, Cedilnik 2018'),
        ('activation_lv_late', 'Activation tardive — Pacing VG',
         'Zones à activation retardée lors d\'un pacing VG (normalisé 0→1)\nSimulation Eikonal, Cedilnik 2018'),
    ]

    def __init__(self, parent=None, previous_config: list = None):
        super().__init__(parent)
        self.setWindowTitle('⊕ Configuration — Score Combiné')
        self.setMinimumWidth(560)
        self.setMinimumHeight(520)
        self._config = {}   # metric_key → weight int 0-100

        # Initialiser depuis config précédente
        if previous_config:
            for entry in previous_config:
                self._config[entry['metric']] = int(entry['weight'] * 100)

        self._build_ui()

    def _build_ui(self):
        root = QtWidgets.QVBoxLayout(self)
        root.setSpacing(8)

        # === Tableau métriques / poids ===
        scroll = QtWidgets.QScrollArea()
        scroll.setWidgetResizable(True)
        scroll.setFrameShape(QtWidgets.QFrame.NoFrame)
        inner = QtWidgets.QWidget()
        grid = QtWidgets.QGridLayout(inner)
        grid.setColumnStretch(1, 1)
        grid.setHorizontalSpacing(10)
        grid.setVerticalSpacing(6)

        # En-têtes
        for col, header in enumerate(['Actif', 'Métrique', 'Poids (%)', '']):
            lbl = QtWidgets.QLabel(f'<b>{header}</b>')
            grid.addWidget(lbl, 0, col)

        self._rows = {}  # metric_key → (checkbox, slider, spinbox)
        for row_idx, (key, label, desc) in enumerate(self.METRIC_ROWS, start=1):
            init_w = self._config.get(key, 0)

            chk = QtWidgets.QCheckBox()
            chk.setChecked(init_w > 0)
            chk.setToolTip(desc)

            lbl_w = QtWidgets.QLabel(label)
            lbl_w.setToolTip(desc)

            slider = QtWidgets.QSlider(QtCore.Qt.Horizontal)
            slider.setRange(0, 100)
            slider.setValue(init_w)
            slider.setEnabled(init_w > 0)

            spin = QtWidgets.QSpinBox()
            spin.setRange(0, 100)
            spin.setValue(init_w)
            spin.setFixedWidth(55)
            spin.setSuffix(' %')
            spin.setEnabled(init_w > 0)

            # Connexions bidirectionnelles slider ↔ spinbox
            slider.valueChanged.connect(spin.setValue)
            spin.valueChanged.connect(slider.setValue)
            # Activer/désactiver slider selon checkbox
            chk.toggled.connect(lambda checked, s=slider, sp=spin, k=key:
                                 self._on_toggle(checked, s, sp, k))

            grid.addWidget(chk,    row_idx, 0, QtCore.Qt.AlignCenter)
            grid.addWidget(lbl_w,  row_idx, 1)
            grid.addWidget(slider, row_idx, 2)
            grid.addWidget(spin,   row_idx, 3)

            self._rows[key] = (chk, slider, spin)

        scroll.setWidget(inner)
        root.addWidget(scroll, 1)

        # === Boutons ===
        btn_row = QtWidgets.QDialogButtonBox(
            QtWidgets.QDialogButtonBox.Ok | QtWidgets.QDialogButtonBox.Cancel
        )
        btn_row.accepted.connect(self.accept)
        btn_row.rejected.connect(self.reject)

        # Bouton reset
        btn_reset = QtWidgets.QPushButton('Réinitialiser')
        btn_reset.clicked.connect(self._reset)
        btn_row.addButton(btn_reset, QtWidgets.QDialogButtonBox.ResetRole)

        root.addWidget(btn_row)

    def _on_toggle(self, checked: bool, slider, spin, key: str):
        slider.setEnabled(checked)
        spin.setEnabled(checked)
        if not checked:
            slider.setValue(0)

    def _reset(self):
        for key, (chk, slider, spin) in self._rows.items():
            chk.setChecked(False)
            slider.setValue(0)
            spin.setValue(0)

    def get_config(self) -> list:
        """Retourne la config sous forme [{'metric': str, 'weight': float}, ...]"""
        result = []
        for key, (chk, slider, spin) in self._rows.items():
            if chk.isChecked() and spin.value() > 0:
                result.append({
                    'metric': key,
                    'weight': spin.value() / 100.0,
                })
        return result
