


def create_dashboard(result_path, plot_path=None):
    
    import os
    import pandas as pd

    results_df=pd.read_csv(result_path)

    import pandas as pd
    import plotly.graph_objects as go
    import plotly.express as px
    import panel as pn

    pn.extension("plotly")


    # ----------------------------
    # Define radar groups
    # ----------------------------
    """ radar_groups = {
    "Annotation-Dependent Metrics": ["ARI", "AMI", "Homogeneity", "Completeness", "V-Measure", "Purity"],
    "Annotation-Independent Spatially Aware Transcriptomic Coherence Vs Standard Transcriptomic Coherence Vs Average Dispersion Metric": ["Silhouette-Spatial","Silhouette","Average-Dispersion"],#"Average-Dispersion"],
    "Annotation-Independent Transcriptomic Coherence vs Spatial Compactness Metrics (I)": ["Silhouette","ASW"],
    "Annotation-Independent Transcriptomic Coherence vs Spatial Compactness Metrics (II)": ["Davies-Bouldin", "CHAOS", "PAS"]
    } """
    radar_groups = {
        "Annotation-Dependent Metrics (Best=Max)": ["ARI", "AMI", "Completeness","Homogeneity", "V-Measure", "Purity"],
        "Annotation-Independent Transcriptomic Coherence and Spatial Compactness Metrics (Best=Max)": ["Silhouette-Spatial","Silhouette", "ASW"],
        "Annotation-Independent Transcriptomic Coherence and Spatial Compactness Metrics (Best=Min)": ["Davies-Bouldin", "CHAOS"],
        "Annotation-Independent Spatial Compactness Metrics and Fragmentation Level (Best=Min)": ["Average-Dispersion","PAS"]
    }
    # Plain-text options for the dropdown
    radar_group_labels = [
        "Annotation-Dependent Metrics (Best=Max)",
        "Annotation-Independent Transcriptomic Coherence and Spatial Compactness Metrics (Best=Max)",
        "Annotation-Independent Transcriptomic Coherence and Spatial Compactness Metrics (Best=Min)",
        "Annotation-Independent Spatial Compactness Metrics and Fragmentation Level (Best=Min)"
    ]
    radar_group_titles = {
        "Annotation-Dependent Metrics (Best=Max)": 
            "Annotation-Dependent Metrics<br>(Best=Max)",

        "Annotation-Independent Transcriptomic Coherence and Spatial Compactness Metrics (Best=Max)": 
            "Annotation-Independent Transcriptomic Coherence and<br>Spatial Compactness Metrics (Best=Max)",

        "Annotation-Independent Transcriptomic Coherence and Spatial Compactness Metrics (Best=Min)": 
            "Annotation-Independent Transcriptomic Coherence and<br>Spatial Compactness Metrics (Best=Min)",

        "Annotation-Independent Spatial Compactness Metrics and Fragmentation Level (Best=Min)": 
            "Annotation-Independent Spatial Compactness Metrics and<br>Fragmentation Level (Best=Min)"
    }


    # ----------------------------
    # Widgets
    # ----------------------------
    model_selector = pn.widgets.MultiChoice(
        name="Select Methods", 
        value=results_df["method"].tolist(),
        options=results_df["method"].tolist()
    )

    radar_group_selector = pn.widgets.Select(
        name="Metric Group",
        #options=list(radar_groups.keys()),
        options=radar_group_labels,
        value="Annotation-Dependent Metrics (Best=Max)",
    )


    #bubble_x = pn.widgets.Select(name="Bubble X Metric", options=list(results_df.columns), value="Silhouette-Spatial")
    #bubble_y = pn.widgets.Select(name="Bubble Y Metric", options=list(results_df.columns), value="Silhouette")
    bubble_x = pn.widgets.Select(name="Scatter X Metric", options=list(results_df.columns), value="Silhouette-Spatial")
    bubble_y = pn.widgets.Select(name="Scatter Y Metric", options=list(results_df.columns), value="Silhouette")

    # ----------------------------
    # Create FigureWidget for radar chart
    # ----------------------------
    fig_radar = go.FigureWidget()

    def update_radar(event=None):
        group_name = radar_group_selector.value
        selected_methods = model_selector.value
        metrics = radar_groups[group_name]
        fig_radar.data = []  # Clear existing traces
        for method in selected_methods:
            values = results_df.loc[results_df["method"] == method, metrics].values.flatten().tolist()
            fig_radar.add_trace(go.Scatterpolar(
                r=values,
                theta=metrics,
                fill='toself',
                name=method
            ))

        # Apply the HTML title mapping
        html_title = radar_group_titles.get(group_name, group_name)   
        fig_radar.update_layout(
            width=700,
            height=700,
            polar=dict(radialaxis=dict(visible=True)),
            showlegend=True,
            legend=dict(
                orientation="h",
                yanchor="bottom",
                y=-0.25,
                xanchor="center",
                x=0.5,
                font=dict(size=11)
            ),
            title=dict(
                text=f"<b>{html_title}</b>",
                font=dict(size=16)
            )#group_name,
        
        )

    # Initial radar chart
    update_radar()

    def safe_hover(event):
        try:
            point = event.values['points'][0]
            curve = point.get('curveNumber', None)
            idx = point.get('pointNumber', None)
            print(f"Hovered curve {curve}, point {idx}")
        except (KeyError, IndexError):
            pass


    # Watch for widget changes
    radar_group_selector.param.watch(update_radar, 'value')
    model_selector.param.watch(update_radar, 'value')

    # ----------------------------
    # Bubble chart and bar chart
    # ----------------------------
    def make_bubble(selected, x_metric, y_metric):
        df_filtered = results_df[results_df["method"].isin(selected)]
        fig = px.scatter(df_filtered, x=x_metric, y=y_metric,
                     #size="exec_time", 
                     color="method",
                     hover_data=["ARI","AMI","Homogeneity","Purity", 
                            "Completeness","V-Measure","Silhouette-Spatial","Average-Dispersion", "Silhouette", "Davies-Bouldin", "CHAOS", "PAS","ASW" ])
    
    
        #fig.update_layout(width=700, height=500,title=f"<b>Trade-off: {x_metric} vs {y_metric} (bubble size = exec_time (s))<b>")
        fig.update_traces(marker=dict(size=14, opacity=0.8,))
        fig.update_layout(width=700, height=500,title=f"<b>Trade-off: {x_metric} vs {y_metric}<b>")


        fig.update_layout(
         #legend_title_text="Method",
            margin=dict(l=40, r=40, t=40, b=40)
        )
        return fig  


    def make_bar(selected):
        df_filtered = results_df[results_df["method"].isin(selected)].copy()
        # Rename columns for clearer legend labels
        df_filtered = df_filtered.rename(columns={
            "exec_time": "Execution Time (s)",
            "peak_memory": "Peak Memory (MB)"
        })

        fig = px.bar(
            df_filtered,
            x="method",
            y=["Execution Time (s)", "Peak Memory (MB)"],
            barmode="group",
            title="<b>Execution Time (in seconds) & Peak Memory (in MB)</b>"
        )

        fig.update_layout(
            width=800,
            height=600,
            legend_title_text="Metric",
            xaxis_title="Method",
            yaxis_title="Value (seconds, MB)",
            font=dict(size=12),
            title_font=dict(size=16)
        )

        return fig

    bubble_panel = pn.bind(make_bubble, model_selector, bubble_x, bubble_y)

    if "exec_time" and "peak_memory" in results_df.columns:
        bar_panel = pn.bind(make_bar, model_selector)
        bar_panel_title=pn.pane.Markdown("## ⏱️ Computational Resource Usage by Method")
    else: 
        bar_panel=None
        bar_panel_title=None



    # ----------------------------
    # Biological interpretation (static part)
    # ----------------------------
    metric_interpretations = {
        "ARI": "Captures faithful domain recovery by measuring agreement between predicted clusters and known tissue domains.",
        "AMI": "Captures the amount of shared information between predicted and reference tissue domains.",
        "Homogeneity": "Relevant for detecting overmixing and consistent reference tissue domain labeling.",
        "Completeness": "Prevents splitting of the same tissue domain.",
        "V-Measure": "Balances avoiding overmixing and undersplitting of tissue domains.",
        "Purity": "Reflects the extent to which each predicted cluster contains spots from only one reference tissue domain.",
        "Silhouette-Spatial": "Relevant for tissues where local microarchitecture matters by balancing spatial adjacency with expression similarity.",
        "Silhouette": "Valuable for datasets where transcriptomic distinction dominates.",
        "Davies-Bouldin": "Relevant for molecularly distinct regions, capturing their low intra-cluster variance and high inter-cluster variance.",
        "CHAOS": "Captures gradual tissue transitions and anatomical boundaries of smooth, contiguous spatial domains.",
        "PAS": "Captures interactions between adjacent cell populations.",
        "ASW": "Beneficial for tissue domains with topologically distinct regions that optimize spatial local density."
    }

    # ----------------------------
    # Function: Build dynamic interpretation table
    # ----------------------------
    import pandas as pd
    import panel as pn

    def make_interpretation_table(selected_methods):
        rows = []
        df_filtered = results_df[results_df["method"].isin(selected_methods)]
        for metric, desc in metric_interpretations.items():
            # Sort by metric (ascending for DBI, else descending)
            ascending = True if metric == "Davies-Bouldin" else False
            top2 = df_filtered.sort_values(metric, ascending=ascending)["method"].head(2).tolist()
            rows.append([metric, ", ".join(top2), desc])

        interp_df = pd.DataFrame(
            rows,
            columns=["Metric Category", "Top 2 Performing Method", "Biological Relevance / Interpretation"]
        )
        return interp_df

    # ----------------------------
    # Panel widget: dynamic, read-only Tabulator
    # ----------------------------
    def make_readonly_interp_table(selected_methods):
        df = make_interpretation_table(selected_methods)
        table = pn.widgets.Tabulator(
            df,
            pagination="remote",
            page_size=13,
            #sizing_mode="stretch_width",
            selectable=False,
            show_index=False
        )

        # Make table fully read-only
        table.configuration = {
            "selectable": False,
            "movableColumns": False,
            "movableRows": False,
            "columnDefaults": {"editable": False},
            "editable": False,
            "clipboard": False,
            "layout": "fitColumns"
        }
        table.disabled = True
        return table

    interp_panel = pn.bind(make_readonly_interp_table, model_selector)

    # ----------------------------
    # 4 Tables of Metric Scores with optional Color Highlighting + proper sorting
    # ----------------------------
    metric_tables = {
        "Annotation-Dependent Metrics": ["ARI", "AMI", "Purity", "Homogeneity", "Completeness", "V-Measure"],
        "Spatially Aware Transcriptomic Metrics": ["Silhouette-Spatial", "Average-Dispersion"],  # or ["Silhouette-Spatial", "ADI"]
        "Spatial Compactness Metrics": ["CHAOS", "PAS", "ASW"],
        "Transcriptomic Coherence Metrics": ["Davies-Bouldin", "Silhouette"]
    }
    def make_metric_table(metrics, selected_methods):
        """
        Create a read-only, fully sortable Tabulator table in Panel 1.8.2
        without HTML/color formatting (numeric columns sort correctly).
        """
        # 1) prepare dataframe and filter by selected methods
        sub_df = results_df[["method"] + [m for m in metrics if m in results_df.columns]].copy()
        if selected_methods is not None:
            sub_df = sub_df[sub_df["method"].isin(selected_methods)].reset_index(drop=True)

        # 2) ensure numeric columns where appropriate
        for m in metrics:
            if m in sub_df.columns:
                sub_df[m] = pd.to_numeric(sub_df[m], errors="coerce")

        # 3) create Tabulator without columns argument
        table = pn.widgets.Tabulator(
            sub_df,
            #sizing_mode="stretch_width",
            pagination="remote",
            page_size=10,
            selectable=False,
            show_index=False
        )

        # 4) assign columns after init
        cols = [{"field": "method", "title": "Method", "headerSort": True}]
        for m in metrics:
            if m not in sub_df.columns:
                continue
            cols.append({
                "field": m,
                "title": m,
                "sorter": "number",
                "hozAlign": "right",
                "editable": False
            })
        table.columns = cols

        # 5) make table fully read-only
        table.configuration = {
            "selectable": False,
            "movableColumns": False,
            "movableRows": False,
            "columnDefaults": {"editable": False},
            "editable": False,
            "clipboard": False,
            "layout": "fitColumns"
        }
        table.disabled = True

        return table



    # Create four tables (as reactive binds so they update when color_toggle or model_selector changes)
    table1 = pn.bind(lambda sel: make_metric_table(metric_tables["Annotation-Dependent Metrics"], sel), model_selector)
    table2 = pn.bind(lambda sel: make_metric_table(metric_tables["Spatial Compactness Metrics"], sel), model_selector)
    table3 = pn.bind(lambda sel: make_metric_table(metric_tables["Transcriptomic Coherence Metrics"], sel), model_selector)
    table4 = pn.bind(lambda sel: make_metric_table(metric_tables["Spatially Aware Transcriptomic Metrics"], sel), model_selector)

    tables_panel = pn.Column(
    
        pn.pane.Markdown("### Annotation-Dependent Metrics"),
        pn.Column(
        pn.pane.Markdown("### Table 1"),
        table1),
        pn.pane.Markdown("### Annotation-Independent Metrics"),
        pn.Row( 
            pn.Column(
            pn.pane.Markdown("### Table 2: Spatial Compactness Metrics"),
        table2),pn.Spacer(width=150),
        pn.Column(
        pn.pane.Markdown("### Table 3: Transcriptomic Coherence Metrics"),
        table3),pn.Spacer(width=150),
        pn.Column(
        pn.pane.Markdown("### Table 4: Spatial Compactness and Transcriptomic Coherence Metrics"),
        table4))
    
    )
    pn.config.raw_css.append("""
    .tabulator-cell div {
        text-align: center;
        font-weight: 500;
        font-size: 13px;
        border-radius: 4px;
    }   
    """)

    # ----------------------------
    # Display tissue cluster images
    # ----------------------------
    import os
    data_name="151673"

    from PIL import Image
    import numpy as np
    import plotly.express as px

    def make_cluster_panel_zoomable(selected_methods, n_cols=5, thumb_size=300, padding=5):
        """
        Display tissue cluster images with interactive zoom.
        Loads images with PIL and converts to NumPy for Plotly.
        """
        panels = []
        method="ground_truth"
        img_path = os.path.join(plot_path, f"{method}.png")
        if os.path.exists(img_path):
            # Load image with PIL and convert to array
                img = np.array(Image.open(img_path))

                # Create Plotly figure
                fig = px.imshow(img)
                fig.update_layout(
                    margin=dict(l=10, r=10, t=30, b=10),
                    width=thumb_size,
                    height=thumb_size,
                )
                fig.update_xaxes(showticklabels=False).update_yaxes(showticklabels=False)
                fig.update_layout(dragmode='zoom')  # enable zoom

                fig_pane = pn.pane.Plotly(fig, sizing_mode="fixed", width=thumb_size, height=thumb_size)

                panels.append(pn.Column(pn.panel(f"### {method}",align="center"), fig_pane, width=thumb_size, margin=(0, padding)))
            
        for method in selected_methods:
            img_path = os.path.join(plot_path, f"{method}.png")
            if os.path.exists(img_path):
            # Load image with PIL and convert to array
                img = np.array(Image.open(img_path))

                # Create Plotly figure
                fig = px.imshow(img)
                fig.update_layout(
                    margin=dict(l=10, r=10, t=30, b=10),
                    width=thumb_size,
                    height=thumb_size,
                )
                fig.update_xaxes(showticklabels=False).update_yaxes(showticklabels=False)
                fig.update_layout(dragmode='zoom')  # enable zoom

                fig_pane = pn.pane.Plotly(fig, sizing_mode="fixed", width=thumb_size, height=thumb_size)

                panels.append(pn.Column(pn.panel(f"### {method}",align="center"), fig_pane, width=thumb_size, margin=(0, padding)))
            else:
                panels.append(pn.Column(pn.panel(f"### {method}",align="center"), pn.pane.Markdown("Image not found"), width=thumb_size, margin=(0, padding)))
    
        # Arrange in rows
        rows = []
        for i in range(0, len(panels), n_cols):
            row = pn.Row(*panels[i:i+n_cols], sizing_mode="stretch_width", margin=(0, padding))
            rows.append(row)

        # Wrap in scrollable column
        return pn.Column(*rows, sizing_mode="stretch_width", margin=(0, padding), scroll=True, max_height=800)
    if plot_path==None:
        cluster_panel=None
        cluster_panel_title=None
    else:
        cluster_panel = pn.bind(make_cluster_panel_zoomable, model_selector)
        cluster_panel=pn.panel(cluster_panel)
        cluster_panel_title=pn.pane.Markdown("## 🖼️ Cluster Visualizations")


    # ----------------------------
    # Final dashboard layout 
    # ----------------------------
    dashboard = pn.Column(
        "# 📊 MultimetricST Dashboard",
        pn.Row(model_selector, margin=(0, 0, 35, 0) ),
        pn.pane.Markdown("## 📈 Visualization of Metric Behaviour"),
        pn.Row( 
        pn.Column(radar_group_selector,
        pn.Row(pn.panel(fig_radar, align="start"),align="start"),margin=(0, 65, 0, 0)),
        pn.Column(
        pn.Row(bubble_x, bubble_y ,margin=(0, 0, 35, 0) ),
        pn.Row(pn.panel(bubble_panel,align="center",margin=(15, 0, 0, 45)), align="end"),margin=(0, 0, 0, 45))
        ,align="center"),
    
        cluster_panel_title,
        cluster_panel,
        pn.pane.Markdown("## 📋 Evaluation Tables of all Metrics"),
        pn.Column("### Biological Interpretation of Metrics",
        interp_panel,),
    
        tables_panel,
    
        bar_panel_title,
        bar_panel,
    
        margin=(35, 35, 35, 35),   # outer padding for full dashboard
        sizing_mode="stretch_both"
    

        )
    return dashboard



