"""
Dashboard Manager - Tableau de bord complet du ventricule
=========================================================
Génère des graphiques détaillés avec matplotlib intégré dans Qt.
Histogrammes, pie charts, radar, box plots, scatter plots, CDF.
Export PDF/PNG disponible.
"""

import numpy as np
from qtpy import QtWidgets, QtCore
from matplotlib.figure import Figure
from matplotlib.gridspec import GridSpec

try:
    from matplotlib.backends.backend_qtagg import FigureCanvasQTAgg as FigureCanvas
except ImportError:
    from matplotlib.backends.backend_qt5agg import FigureCanvasQTAgg as FigureCanvas


# ═══════════════════════════════════════════════════════════════════════════════
# Couleurs médicales cohérentes
# ═══════════════════════════════════════════════════════════════════════════════
COLORS = {
    'scar':       '#c62828',
    'border':     '#ef6c00',
    'healthy':    '#2e7d32',
    'accent':     '#1565c0',
    'title':      '#1a237e',
    'text':       '#333333',
    'light_text': '#777777',
    'bg_card':    '#f5f5f5',
    'grid':       '#e0e0e0',
}


class DashboardWindow(QtWidgets.QDialog):
    """Fenêtre de tableau de bord complet avec graphiques matplotlib"""

    def __init__(self, data_manager, parent=None):
        super().__init__(parent)
        self.dm = data_manager
        self.setWindowTitle("LV Explorer \u2014 Ventricular Dashboard Report")
        self.setMinimumSize(1400, 900)
        self.resize(1600, 1050)

        # Pré-calculer toutes les métriques
        self.dm.compute_all_metrics()
        self.stats = self.dm.get_statistics()

        self._setup_ui()

    # ─────────────────────────────────────────────────────────────────────────
    # UI
    # ─────────────────────────────────────────────────────────────────────────
    def _setup_ui(self):
        layout = QtWidgets.QVBoxLayout(self)

        # Toolbar
        toolbar = QtWidgets.QHBoxLayout()

        export_pdf = QtWidgets.QPushButton("Export PDF")
        export_pdf.setStyleSheet(
            "QPushButton{background:#d32f2f;color:white;padding:8px 20px;"
            "border-radius:4px;font-weight:bold}"
            "QPushButton:hover{background:#b71c1c}")
        export_pdf.clicked.connect(self._export_pdf)
        toolbar.addWidget(export_pdf)

        export_png = QtWidgets.QPushButton("Export PNG")
        export_png.setStyleSheet(
            "QPushButton{background:#1565c0;color:white;padding:8px 20px;"
            "border-radius:4px;font-weight:bold}"
            "QPushButton:hover{background:#0d47a1}")
        export_png.clicked.connect(self._export_png)
        toolbar.addWidget(export_png)

        toolbar.addStretch()

        close_btn = QtWidgets.QPushButton("Close")
        close_btn.clicked.connect(self.close)
        toolbar.addWidget(close_btn)

        layout.addLayout(toolbar)

        # Scroll area contenant la figure
        scroll = QtWidgets.QScrollArea()
        scroll.setWidgetResizable(True)

        self.figure = self._create_dashboard_figure()
        self.canvas = FigureCanvas(self.figure)
        self.canvas.setMinimumHeight(2600)

        scroll.setWidget(self.canvas)
        layout.addWidget(scroll)

    # ─────────────────────────────────────────────────────────────────────────
    # Génération de la figure principale
    # ─────────────────────────────────────────────────────────────────────────
    def _create_dashboard_figure(self) -> Figure:
        fig = Figure(figsize=(18, 28), dpi=100, facecolor='white')
        fig.suptitle('LV Explorer \u2014 Ventricular Dashboard Report',
                     fontsize=22, fontweight='bold', y=0.997, color=COLORS['title'])

        gs = GridSpec(5, 3, figure=fig, hspace=0.45, wspace=0.35,
                      top=0.975, bottom=0.025, left=0.06, right=0.97)

        # Row 0 : Summary + Composition
        self._plot_summary_card(fig, gs[0, 0])
        self._plot_tissue_pie(fig, gs[0, 1])
        self._plot_scar_composition(fig, gs[0, 2])

        # Row 1 : Wall Thickness
        self._plot_thickness_histogram(fig, gs[1, 0])
        self._plot_thickness_boxplot(fig, gs[1, 1])
        self._plot_thickness_cumulative(fig, gs[1, 2])

        # Row 2 : Arrhythmia Substrate
        self._plot_ciaccio_histogram(fig, gs[2, 0])
        self._plot_entropy_histogram(fig, gs[2, 1])
        self._plot_channels_histogram(fig, gs[2, 2])

        # Row 3 : Scar & Transmurality
        self._plot_transmurality_histogram(fig, gs[3, 0])
        self._plot_transmurality_pie(fig, gs[3, 1])
        self._plot_risk_radar(fig, gs[3, 2])

        # Row 4 : Additional analyses
        self._plot_regional_bars(fig, gs[4, 0])
        self._plot_scar_characteristics(fig, gs[4, 1])
        self._plot_detailed_summary(fig, gs[4, 2])

        return fig

    # ═══════════════════════════════════════════════════════════════════════════
    # ROW 0 : Summary
    # ═══════════════════════════════════════════════════════════════════════════
    def _plot_summary_card(self, fig, gs):
        ax = fig.add_subplot(gs)
        ax.axis('off')
        ax.set_title('Key Metrics', fontsize=14, fontweight='bold',
                     color=COLORS['title'], loc='left')

        s = self.stats
        lines = []

        # Geometry
        lines.append(('GEOMETRY', '', True))
        if 'LV_Surface_cm2' in s:
            lines.append(('  Surface', f"{s['LV_Surface_cm2']:.1f} cm\u00b2", False))
        if 'LV_Volume_mL' in s:
            lines.append(('  Volume', f"{s['LV_Volume_mL']:.1f} mL", False))
        if 'LV_Long_Axis_mm' in s:
            lines.append(('  Long Axis', f"{s['LV_Long_Axis_mm']:.1f} mm", False))
        if 'LV_Sphericity' in s:
            lines.append(('  Sphericity', f"{s['LV_Sphericity']:.3f}", False))
        if 'LV_N_Points' in s:
            lines.append(('  Mesh Points', f"{s['LV_N_Points']:,}", False))

        lines.append(('', '', False))
        lines.append(('THICKNESS', '', True))
        if 'T_mean' in s:
            lines.append(('  Mean \u00b1 SD',
                          f"{s['T_mean']:.2f} \u00b1 {s['T_std']:.2f} mm", False))
            lines.append(('  Range',
                          f"[{s['T_min']:.2f} \u2013 {s['T_max']:.2f}] mm", False))
            lines.append(('  Median [IQR]',
                          f"{s['T_median']:.2f} [{s.get('T_IQR', 0):.2f}] mm", False))

        lines.append(('', '', False))
        lines.append(('SCAR', '', True))
        if 'Dense_Scar_Surface_cm2' in s:
            lines.append(('  Dense Scar',
                          f"{s['Dense_Scar_Surface_cm2']:.2f} cm\u00b2", False))
            if 'Dense_Scar_Burden_Pct' in s:
                lines.append(('  Scar Burden',
                              f"{s['Dense_Scar_Burden_Pct']:.1f}%", False))
        if 'Dense_Scar_Tissue_Volume_cm3' in s:
            lines.append(('  Tissue Vol.',
                          f"{s['Dense_Scar_Tissue_Volume_cm3']:.2f} cm\u00b3", False))

        y = 0.96
        for label, value, is_header in lines:
            if is_header:
                ax.text(0.02, y, label, transform=ax.transAxes,
                        fontsize=10, fontweight='bold', color=COLORS['title'],
                        fontfamily='monospace')
            elif label:
                ax.text(0.04, y, label, transform=ax.transAxes,
                        fontsize=8.5, color=COLORS['light_text'],
                        fontfamily='monospace')
                ax.text(0.95, y, value, transform=ax.transAxes,
                        fontsize=8.5, fontweight='bold', color=COLORS['text'],
                        fontfamily='monospace', ha='right')
            y -= 0.050

    def _plot_tissue_pie(self, fig, gs):
        ax = fig.add_subplot(gs)
        ax.set_title('Tissue Composition',
                     fontsize=14, fontweight='bold', color=COLORS['title'])

        s = self.stats
        if 'Scar_Pct' in s:
            sizes = [s['Scar_Pct'], s['Border_Zone_Pct'], s['Healthy_Pct']]
            labels = [f"Scar (<1mm)\n{sizes[0]:.1f}%",
                      f"Border (1-5mm)\n{sizes[1]:.1f}%",
                      f"Healthy (>5mm)\n{sizes[2]:.1f}%"]
            colors = [COLORS['scar'], COLORS['border'], COLORS['healthy']]
            explode = (0.05, 0.02, 0)

            ax.pie(sizes, labels=labels, colors=colors, explode=explode,
                   startangle=90, textprops={'fontsize': 9},
                   wedgeprops={'edgecolor': 'white', 'linewidth': 1.5})
            ax.axis('equal')
        else:
            self._no_data(ax, 'No thickness data')

    def _plot_scar_composition(self, fig, gs):
        ax = fig.add_subplot(gs)
        ax.set_title('Scar Analysis (à revoir)',
                     fontsize=14, fontweight='bold', color=COLORS['title'])

        s = self.stats
        cats, vals, cols = [], [], []

        mapping = [
            ('Dense_Scar_Surface_cm2', 'Dense\nScar',  '#c62828'),
            ('Scar_LE_Surface_cm2',    'Scar\n(LE)',   '#e53935'),
            ('Scar_Endo_Surface_cm2',  'Endo',         '#ef5350'),
            ('Scar_Intra_Surface_cm2', 'Intra',        '#ff7043'),
            ('Scar_Epi_Surface_cm2',   'Epi',          '#ffa726'),
            ('Fat_Surface_cm2',        'Fat',          '#ffd54f'),
        ]
        for key, label, color in mapping:
            if key in s:
                cats.append(label); vals.append(s[key]); cols.append(color)

        if cats:
            bars = ax.bar(cats, vals, color=cols, edgecolor='white', linewidth=0.5)
            ax.set_ylabel('Surface (cm\u00b2)', fontsize=10)
            self._clean_axes(ax)
            for bar, val in zip(bars, vals):
                ax.text(bar.get_x() + bar.get_width() / 2,
                        bar.get_height() + max(vals) * 0.02,
                        f'{val:.2f}', ha='center', va='bottom',
                        fontsize=8, fontweight='bold')
        else:
            self._no_data(ax, 'No scar data')

    # ═══════════════════════════════════════════════════════════════════════════
    # ROW 1 : Wall Thickness
    # ═══════════════════════════════════════════════════════════════════════════
    def _plot_thickness_histogram(self, fig, gs):
        ax = fig.add_subplot(gs)
        ax.set_title('Wall Thickness Distribution',
                     fontsize=14, fontweight='bold', color=COLORS['title'])

        T = self.dm.get_scalar_data('wall_thickness', 'Wall_Thickness')
        if T is None:
            self._no_data(ax); return

        ax.hist(T, bins=60, color=COLORS['accent'], alpha=0.7,
                edgecolor='white', linewidth=0.4, density=True)

        # KDE
        try:
            from scipy.stats import gaussian_kde
            kde = gaussian_kde(T, bw_method=0.15)
            x = np.linspace(max(T.min(), 0), T.max(), 300)
            ax.plot(x, kde(x), color='#0d47a1', linewidth=2, label='KDE')
        except Exception:
            pass

        # Seuils
        ax.axvline(1.0, color=COLORS['scar'], ls='--', lw=1.5,
                   label='Scar (1 mm)')
        ax.axvline(5.0, color=COLORS['border'], ls='--', lw=1.5,
                   label='Border (5 mm)')
        mean_v = float(np.mean(T))
        ax.axvline(mean_v, color=COLORS['healthy'], ls='-', lw=1.5,
                   label=f'Mean ({mean_v:.2f} mm)')

        ax.set_xlabel('Thickness (mm)', fontsize=10)
        ax.set_ylabel('Density', fontsize=10)
        ax.legend(fontsize=8, loc='upper right')
        self._clean_axes(ax)

    def _plot_thickness_boxplot(self, fig, gs):
        ax = fig.add_subplot(gs)
        ax.set_title('Thickness Statistics',
                     fontsize=14, fontweight='bold', color=COLORS['title'])

        T = self.dm.get_scalar_data('wall_thickness', 'Wall_Thickness')
        if T is None:
            self._no_data(ax); return

        bp = ax.boxplot(
            T, vert=True, patch_artist=True,
            boxprops=dict(facecolor='#bbdefb', linewidth=1.5),
            medianprops=dict(color=COLORS['scar'], linewidth=2),
            whiskerprops=dict(linewidth=1.5),
            capprops=dict(linewidth=1.5),
            flierprops=dict(marker='o', markerfacecolor='#e53935',
                            markersize=3, alpha=0.3))

        ax.set_ylabel('Thickness (mm)', fontsize=10)
        ax.set_xticklabels(['Wall\nThickness'])

        s = self.stats
        if 'T_mean' in s:
            ax.axhline(s['T_mean'], color=COLORS['healthy'], ls='--', alpha=0.5)
            offset = 1.2
            annotations = [
                (s['T_mean'], f"Mean: {s['T_mean']:.2f}", COLORS['healthy']),
                (s['T_median'], f"Median: {s['T_median']:.2f}", COLORS['scar']),
                (s.get('T_p5', 0), f"P5: {s.get('T_p5', 0):.2f}", '#999'),
                (s.get('T_p95', 0), f"P95: {s.get('T_p95', 0):.2f}", '#999'),
            ]
            for yval, txt, col in annotations:
                ax.text(offset, yval, txt, fontsize=8, color=col, va='center')
        self._clean_axes(ax)

    def _plot_thickness_cumulative(self, fig, gs):
        ax = fig.add_subplot(gs)
        ax.set_title('Cumulative Distribution (CDF)',
                     fontsize=14, fontweight='bold', color=COLORS['title'])

        T = self.dm.get_scalar_data('wall_thickness', 'Wall_Thickness')
        if T is None:
            self._no_data(ax); return

        sorted_T = np.sort(T)
        cdf = np.arange(1, len(sorted_T) + 1) / len(sorted_T) * 100

        ax.plot(sorted_T, cdf, color=COLORS['accent'], linewidth=2)
        ax.fill_between(sorted_T, cdf, alpha=0.08, color=COLORS['accent'])

        ax.axvline(1.0, color=COLORS['scar'], ls='--', alpha=0.5)
        ax.axvline(5.0, color=COLORS['border'], ls='--', alpha=0.5)

        for p, col in [(5, '#aaa'), (25, '#888'), (50, '#555'),
                        (75, '#888'), (95, '#aaa')]:
            val = float(np.percentile(T, p))
            ax.axhline(p, color=col, ls=':', alpha=0.25)
            ax.plot(val, p, 'o', color=COLORS['scar'], markersize=5)
            ax.text(val + 0.15, p + 2.5, f'P{p}={val:.1f}',
                    fontsize=7, color=COLORS['text'])

        ax.set_xlabel('Thickness (mm)', fontsize=10)
        ax.set_ylabel('Cumulative %', fontsize=10)
        ax.set_ylim(0, 105)
        self._clean_axes(ax)
        ax.grid(True, alpha=0.15)

    # ═══════════════════════════════════════════════════════════════════════════
    # ROW 2 : Arrhythmia Substrate
    # ═══════════════════════════════════════════════════════════════════════════
    def _plot_ciaccio_histogram(self, fig, gs):
        ax = fig.add_subplot(gs)
        ax.set_title('Ciaccio Ratio (\u03c1) Distribution',
                     fontsize=14, fontweight='bold', color=COLORS['title'])

        rho = self.dm.get_scalar_data('ciaccio_ratio', 'Ciaccio_Ratio')
        if rho is None:
            self._no_data(ax); return

        rho_disp = rho[rho < 2.0]
        ax.hist(rho_disp, bins=60, color='#ff7043', alpha=0.7,
                edgecolor='white', linewidth=0.4, density=True)

        s = self.stats
        if 'Rho_p95' in s:
            ax.axvline(s['Rho_p95'], color='red', ls='--', lw=2,
                       label=f"95th pct ({s['Rho_p95']:.3f})")
        if 'Rho_Mean' in s:
            ax.axvline(s['Rho_Mean'], color=COLORS['healthy'], ls='-', lw=1.5,
                       label=f"Mean ({s['Rho_Mean']:.3f})")

        ax.set_xlabel('Ciaccio Ratio (\u03c1)', fontsize=10)
        ax.set_ylabel('Density', fontsize=10)
        ax.legend(fontsize=8)
        self._clean_axes(ax)

    def _plot_channels_histogram(self, fig, gs):
        ax = fig.add_subplot(gs)
        ax.set_title('Channel Score Distribution',
                     fontsize=14, fontweight='bold', color=COLORS['title'])

        ch_score = self.dm.get_scalar_data('channels', 'Channel_Score')
        if ch_score is None:
            self._no_data(ax, 'Not computed'); return

        ch_active = ch_score[ch_score > 0]
        if len(ch_active) == 0:
            self._no_data(ax, 'No channels detected'); return

        ax.hist(ch_active, bins=40, color='#1976d2', alpha=0.7,
                edgecolor='white', linewidth=0.4, density=True)
        m = float(np.mean(ch_active))
        ax.axvline(m, color=COLORS['healthy'], ls='-', alpha=0.7,
                   label=f'Mean ({m:.3f})')

        ax.set_xlabel('Channel Score', fontsize=10)
        ax.set_ylabel('Density', fontsize=10)
        ax.legend(fontsize=8)
        self._clean_axes(ax)

    def _plot_entropy_histogram(self, fig, gs):
        ax = fig.add_subplot(gs)
        ax.set_title('Local Entropy Distribution',
                     fontsize=14, fontweight='bold', color=COLORS['title'])

        ent = self.dm.get_scalar_data('local_entropy', 'Local_Entropy')
        if ent is None:
            self._no_data(ax, 'Not computed'); return

        ax.hist(ent, bins=40, color='#5c6bc0', alpha=0.7,
                edgecolor='white', linewidth=0.4, density=True)
        m = float(np.mean(ent))
        p75 = float(np.percentile(ent, 75))
        p95 = float(np.percentile(ent, 95))
        ax.axvline(m, color=COLORS['healthy'], ls='-', alpha=0.7,
                   label=f'Mean ({m:.3f})')
        ax.axvline(p75, color='#ff9800', ls='--', alpha=0.6,
                   label=f'P75 ({p75:.3f})')
        ax.axvline(p95, color='#c62828', ls='--', alpha=0.6,
                   label=f'P95 ({p95:.3f})')

        ax.set_xlabel('Entropy (normalized)', fontsize=10)
        ax.set_ylabel('Density', fontsize=10)
        ax.legend(fontsize=8)
        self._clean_axes(ax)

    # ═══════════════════════════════════════════════════════════════════════════
    # ROW 3 : Scar & Risk
    # ═══════════════════════════════════════════════════════════════════════════
    def _plot_transmurality_histogram(self, fig, gs):
        ax = fig.add_subplot(gs)
        ax.set_title('Scar Transmurality',
                     fontsize=14, fontweight='bold', color=COLORS['title'])

        trans = self.dm.get_scalar_data('transmurality', 'Transmurality')
        if trans is None:
            self._no_data(ax, 'No transmurality data'); return

        ax.hist(trans, bins=40, color='#ef5350', alpha=0.7,
                edgecolor='white', linewidth=0.4)
        for val, label, col in [(25, 'Sub-endo/Mid', '#ff9800'),
                                 (50, 'Mid/Sub-epi', '#1976d2'),
                                 (75, 'Sub-epi/Trans', '#7b1fa2')]:
            ax.axvline(val, color=col, ls='--', alpha=0.7, label=f'{label} ({val}%)')

        ax.set_xlabel('Transmurality (%)', fontsize=10)
        ax.set_ylabel('Count', fontsize=10)
        ax.legend(fontsize=7)
        self._clean_axes(ax)

    def _plot_transmurality_pie(self, fig, gs):
        ax = fig.add_subplot(gs)
        ax.set_title('Transmurality Classification',
                     fontsize=14, fontweight='bold', color=COLORS['title'])

        s = self.stats
        keys = ['Trans_None_Pct', 'Trans_Subendo_Pct', 'Trans_Midmural_Pct',
                'Trans_Subepi_Pct', 'Trans_Transmural_Pct']

        if all(k in s for k in keys):
            sizes  = [s[k] for k in keys]
            labels = ['None', 'Sub-endo\n(0-25%)', 'Midmural\n(25-50%)',
                      'Sub-epi\n(50-75%)', 'Transmural\n(>75%)']
            colors = ['#e8f5e9', '#fff9c4', '#ffe0b2', '#ffccbc', '#ef9a9a']

            non_zero = [(sz, lb, cl) for sz, lb, cl
                        in zip(sizes, labels, colors) if sz > 0.5]
            if non_zero:
                sz, lb, cl = zip(*non_zero)
                lb_fmt = [f'{l}\n{v:.1f}%' for l, v in zip(lb, sz)]
                ax.pie(sz, labels=lb_fmt, colors=cl, startangle=90,
                       textprops={'fontsize': 8},
                       wedgeprops={'edgecolor': 'white', 'linewidth': 1})
                ax.axis('equal')
            else:
                self._no_data(ax, 'No transmural scar')
        else:
            self._no_data(ax, 'No transmurality data')

    def _plot_risk_radar(self, fig, gs):
        """Radar chart du profil de risque multi-axes."""
        ax = fig.add_subplot(gs, polar=True)
        ax.set_title('Risk Profile',
                     fontsize=14, fontweight='bold', color=COLORS['title'], pad=25)

        s = self.stats
        risk = s.get('_risk', {})
        if not risk:
            ax.set_rticks([]); ax.set_xticks([])
            ax.text(0.5, 0.5, 'Insufficient data\nfor risk profile',
                    ha='center', va='center', fontsize=10, color='gray',
                    transform=ax.transAxes)
            return

        # Axes du radar (sans mentions d'articles)
        radar_axes = [
            ('Scar\nBurden',    'scar_burden_norm'),
            ('Scar\nVolume',    'scar_volume_norm'),
            ('Decel.\nZones',   'dz_extent_norm'),
            ('Wall\nThinning',  'thinning_norm'),
            ('Trans-\nmurality','transmural_norm'),
            ('Sphericity',      'sphericity_norm'),
            ('Channels',        'channel_extent_norm'),
            ('Isthmus',         'isthmus_norm'),
            ('Entropy',         'entropy_norm'),
        ]

        metrics, values = [], []
        for label, key in radar_axes:
            if key in risk:
                metrics.append(label)
                values.append(risk[key])

        if len(metrics) >= 3:
            N = len(metrics)
            angles = np.linspace(0, 2 * np.pi, N, endpoint=False).tolist()
            v_plot = values + [values[0]]
            angles_plot = angles + [angles[0]]

            ax.plot(angles_plot, v_plot, 'o-', linewidth=2,
                    color='#e53935', markersize=6)
            ax.fill(angles_plot, v_plot, alpha=0.20, color='#e53935')

            ax.set_xticks(angles)
            ax.set_xticklabels(metrics, fontsize=7)
            ax.set_ylim(0, 1.05)
            ax.set_yticks([0.25, 0.5, 0.75, 1.0])
            ax.set_yticklabels(['Low', 'Mod', 'High', 'V.High'],
                               fontsize=7, color='gray')
            ax.grid(True, alpha=0.3)

            # Score composite en bas
            score = risk.get('composite_score', 0)
            level = risk.get('composite_level', '?')
            col = ('#2e7d32' if score < 0.25 else '#ef6c00'
                   if score < 0.50 else '#c62828' if score < 0.75 else '#4a148c')
            ax.text(0.5, -0.12, f'Composite: {score:.2f} ({level})',
                    ha='center', fontsize=10, fontweight='bold', color=col,
                    transform=ax.transAxes)
        else:
            ax.set_rticks([])
            ax.set_xticks([])
            ax.text(0.5, 0.5, 'Insufficient data\nfor risk profile',
                    ha='center', va='center', fontsize=10, color='gray',
                    transform=ax.transAxes)

    # ═══════════════════════════════════════════════════════════════════════════
    # ROW 4 : Additional
    # ═══════════════════════════════════════════════════════════════════════════
    def _plot_regional_bars(self, fig, gs):
        ax = fig.add_subplot(gs)
        ax.set_title('Surface Area Breakdown',
                     fontsize=14, fontweight='bold', color=COLORS['title'])

        s = self.stats
        cats, vals, cols = [], [], []

        if 'Scar_Area_cm2' in s:
            for key, label, color in [
                ('Scar_Area_cm2',    'Scar\n(<1mm)',    COLORS['scar']),
                ('Border_Area_cm2',  'Border\n(1-5mm)', COLORS['border']),
                ('Healthy_Area_cm2', 'Healthy\n(>5mm)', COLORS['healthy']),
            ]:
                if key in s:
                    cats.append(label); vals.append(s[key]); cols.append(color)

        if cats:
            bars = ax.bar(cats, vals, color=cols, edgecolor='white', linewidth=0.5)
            ax.set_ylabel('Surface (cm\u00b2)', fontsize=10)
            self._clean_axes(ax)
            for bar, val in zip(bars, vals):
                ax.text(bar.get_x() + bar.get_width() / 2,
                        bar.get_height() + max(vals) * 0.02,
                        f'{val:.1f}', ha='center', va='bottom',
                        fontsize=9, fontweight='bold')
        else:
            self._no_data(ax, 'No area data')

    def _plot_scar_characteristics(self, fig, gs):
        """Caractéristiques du substrat."""
        ax = fig.add_subplot(gs)
        ax.set_title('Substrate Characteristics',
                     fontsize=14, fontweight='bold', color=COLORS['title'])
        ax.axis('off')

        s = self.stats
        lines = []

        lines.append(('SCAR MORPHOLOGY', COLORS['title'], True))
        if 'Dense_Scar_Surface_cm2' in s:
            lines.append((f"  Surface:     {s['Dense_Scar_Surface_cm2']:>7.2f} cm2",
                           COLORS['text'], False))
        if 'Dense_Scar_Burden_Pct' in s:
            v = s['Dense_Scar_Burden_Pct']
            c = '#c62828' if v >= 10 else COLORS['text']
            flag = ' (!)' if v >= 10 else ''
            lines.append((f"  Burden:      {v:>7.1f}%{flag}", c, False))
        if 'Dense_Scar_Tissue_Volume_cm3' in s:
            v = s['Dense_Scar_Tissue_Volume_cm3']
            c = '#c62828' if v >= 37.3 else COLORS['text']
            flag = ' (!)' if v >= 37.3 else ''
            lines.append((f"  Volume:      {v:>7.2f} cm3{flag}", c, False))
        if 'Scar_Mean_Thickness_mm' in s:
            lines.append((f"  Mean T scar: {s['Scar_Mean_Thickness_mm']:>7.2f} mm",
                           COLORS['text'], False))
        if 'Scar_Compactness' in s:
            lines.append((f"  Compactness: {s['Scar_Compactness']:>7.3f}",
                           COLORS['text'], False))

        lines.append(('', COLORS['text'], False))
        lines.append(('ARRHYTHMIA SUBSTRATE', COLORS['title'], True))
        if 'DZ_Area_Pct' in s:
            lines.append((f"  DZ extent:   {s['DZ_Area_Pct']:>7.1f}%",
                           COLORS['text'], False))
        if 'Channel_Pct' in s:
            lines.append((f"  Channel ext: {s['Channel_Pct']:>7.1f}%",
                           COLORS['text'], False))
        if 'Channel_Score_Mean' in s:
            lines.append((f"  Ch. score:   {s['Channel_Score_Mean']:>7.3f}",
                           COLORS['text'], False))
        if 'Isthmus_N_Points' in s:
            lines.append((f"  Isthmus pts: {s['Isthmus_N_Points']:>7d}",
                           COLORS['text'], False))

        lines.append(('', COLORS['text'], False))
        lines.append(('HETEROGENEITY', COLORS['title'], True))
        if 'Entropy_Mean' in s:
            lines.append((f"  Entropy avg: {s['Entropy_Mean']:>7.4f}",
                           COLORS['text'], False))
        if 'Entropy_P75' in s:
            lines.append((f"  Entropy P75: {s['Entropy_P75']:>7.4f}",
                           COLORS['text'], False))
        if 'Entropy_P95' in s:
            lines.append((f"  Entropy P95: {s['Entropy_P95']:>7.4f}",
                           COLORS['text'], False))

        y = 0.95
        for text, color, bold in lines:
            ax.text(0.04, y, text, transform=ax.transAxes,
                    fontsize=7.5, fontfamily='monospace', color=color,
                    fontweight='bold' if bold else 'normal')
            y -= 0.050

    def _plot_detailed_summary(self, fig, gs):
        """Risk Assessment sans citations littérature."""
        ax = fig.add_subplot(gs)
        ax.axis('off')
        ax.set_title('Risk Assessment',
                     fontsize=14, fontweight='bold', color=COLORS['title'],
                     loc='left')

        s = self.stats
        risk = s.get('_risk', {})

        lines = []  # (text, color, bold)

        # --- Threshold-based metrics ---
        lines.append(('THRESHOLD-BASED METRICS', COLORS['title'], True))

        sb = risk.get('scar_burden_pct', 0)
        sb_flag = '!' if sb >= 10 else ' '
        lines.append((f"  {sb_flag} Scar burden:    {sb:>6.1f}%  (>= 10%)",
                       '#c62828' if sb >= 10 else '#2e7d32', False))

        sv = risk.get('scar_volume_cm3', 0)
        sv_flag = '!' if sv >= 37.3 else ' '
        lines.append((f"  {sv_flag} Scar volume:   {sv:>6.1f} cm3 (>= 37.3)",
                       '#c62828' if sv >= 37.3 else '#2e7d32', False))

        dz = risk.get('dz_extent_pct', 0)
        lines.append((f"    DZ extent:     {dz:>6.1f}%",
                       COLORS['text'], False))

        sph = risk.get('sphericity', 0)
        sph_flag = '!' if sph >= 0.70 else ' '
        lines.append((f"  {sph_flag} Sphericity:    {sph:>6.3f}  (>= 0.70)",
                       '#c62828' if sph >= 0.70 else '#2e7d32', False))

        lines.append(('', COLORS['text'], False))

        # --- Continuous metrics ---
        lines.append(('CONTINUOUS METRICS',
                       COLORS['title'], True))

        ch = risk.get('channel_extent_pct', 0)
        lines.append((f"    Channel extent:  {ch:>6.1f}%", COLORS['text'], False))

        ip = risk.get('isthmus_pct', 0)
        lines.append((f"    Isthmus:         {ip:>6.1f}%", COLORS['text'], False))

        em = risk.get('entropy_mean', 0)
        lines.append((f"    Entropy mean:    {em:>6.4f}", COLORS['text'], False))

        tp = risk.get('transmural_pct', 0)
        lines.append((f"    Transmural scar: {tp:>6.1f}%", COLORS['text'], False))

        thin = risk.get('thinning_pct', 0)
        lines.append((f"    Wall thinning:   {thin:>6.1f}%", COLORS['text'], False))

        lines.append(('', COLORS['text'], False))

        # --- Composite score ---
        score = risk.get('composite_score', 0)
        level = risk.get('composite_level', '?')
        col = ('#2e7d32' if score < 0.25 else '#ef6c00'
               if score < 0.50 else '#c62828' if score < 0.75 else '#4a148c')
        lines.append(('COMPOSITE RISK SCORE', COLORS['title'], True))
        lines.append((f"  Score: {score:.3f} / 1.000 -> {level}", col, True))
        lines.append(('  EXPLORATORY - not clinically validated',
                       '#c62828', False))

        y = 0.95
        for text, color, bold in lines:
            ax.text(0.03, y, text, transform=ax.transAxes,
                    fontsize=7, fontfamily='monospace', color=color,
                    fontweight='bold' if bold else 'normal')
            y -= 0.052

    # ─────────────────────────────────────────────────────────────────────────
    # Helpers
    # ─────────────────────────────────────────────────────────────────────────
    @staticmethod
    def _clean_axes(ax):
        ax.spines['top'].set_visible(False)
        ax.spines['right'].set_visible(False)

    @staticmethod
    def _no_data(ax, msg='No data'):
        ax.text(0.5, 0.5, msg, ha='center', va='center',
                fontsize=12, color='gray', transform=ax.transAxes)
        ax.axis('off')

    def _export_pdf(self):
        path, _ = QtWidgets.QFileDialog.getSaveFileName(
            self, "Export Dashboard PDF", "LV_Dashboard_Report.pdf",
            "PDF Files (*.pdf)")
        if path:
            self.figure.savefig(path, format='pdf', dpi=150, bbox_inches='tight')
            QtWidgets.QMessageBox.information(
                self, "Export", f"Dashboard exported to:\n{path}")

    def _export_png(self):
        path, _ = QtWidgets.QFileDialog.getSaveFileName(
            self, "Export Dashboard PNG", "LV_Dashboard_Report.png",
            "PNG Files (*.png)")
        if path:
            self.figure.savefig(path, format='png', dpi=200, bbox_inches='tight')
            QtWidgets.QMessageBox.information(
                self, "Export", f"Dashboard exported to:\n{path}")
