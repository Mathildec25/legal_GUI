from dash import html, dash_table
from utils.data_loader import df

# Function that generates the "Database" page
def get_database_page():
    return html.Div([

        # Page title
        html.H2("Full Database", style={
            "marginTop": "20px",
            "marginBottom": "20px",
            "marginRight": "100px",
            "textAlign": "center",
            "fontSize": "36px",
            "fontWeight": "bold",
            "color": "#003366",
            "fontFamily": "Arial, sans-serif"
        }),

        # Interactive data table
        html.Div([
            dash_table.DataTable(
                id="database-table",
                columns=[{"name": col, "id": col} for col in df.columns],
                data=df.to_dict("records"),

                # Table is not editable by the user
                editable=False,

                # Native filtering and sorting
                filter_action="native",
                sort_action="native",

                # Pagination disabled for smooth scrolling
                page_action="none",
                fixed_rows={"headers": True},  # Sticky header on scroll

                # Global style of the table container
                style_table={
                    "width": "calc(100% - 80px)",       # Adjusted to fit next to the sidebar
                    "height": "calc(100vh - 160px)",    # Full visible height
                    "overflowX": "auto",
                    "overflowY": "auto",
                    "marginLeft": "0px",
                    "marginRight": "20px",
                    "border": "1px solid #ddd",
                    "boxShadow": "0 1px 3px rgba(0, 0, 0, 0.1)"
                },

                # Style for all cells
                style_cell={
                    "textAlign": "left",
                    "padding": "10px 14px",
                    "minWidth": "120px",
                    "maxWidth": "300px",
                    "whiteSpace": "normal",
                    "border": "1px solid #f0f0f0",
                    "overflow": "hidden",
                    "textOverflow": "ellipsis",
                    "fontFamily": "Segoe UI, sans-serif",
                    "fontSize": "14px"
                },

                # Header style (sticky + yellow)
                style_header={
                    "backgroundColor": "#ffcb47",
                    "color": "black",
                    "fontWeight": "bold",
                    "position": "sticky",
                    "top": 0,
                    "zIndex": 1,
                    "border": "1px solid #ddd"
                },

                # Data style
                style_data={
                    "backgroundColor": "white",
                    "color": "black"
                },

                # Highlighting: odd rows and active row
                style_data_conditional=[
                    {
                        "if": {"row_index": "odd"},
                        "backgroundColor": "#fafafa"
                    },
                    {
                        "if": {"state": "active"},
                        "backgroundColor": "#d7f0e2",
                        "border": "1px solid #999"
                    }
                ],

                # Performance + CSV export
                virtualization=True,
                export_format="csv"
            )
        ])
    ],
    style={
        "width": "100%",
        "overflow": "hidden",
        "height": "100vh"
    })
