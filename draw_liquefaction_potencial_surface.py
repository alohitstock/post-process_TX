import pathlib
import pandas as pd
import numpy as np
import plotly.graph_objects as go

def main():
    csv_path = r"D:\OneDrive\01_Research\01_Doctor\07_Triaxial_Test\07_Result\2020\20200122\20200122_15_liq_CSR015-OCR1-sc60\20200122_15_liq_CSR015-OCR1-sc60_postprocessed.csv"
    csv_path = pathlib.Path(csv_path)
    csv_dir = csv_path.parent

    data_csv = pd.read_csv(csv_path)

    # print(data_csv.columns)
    # data_x = 1 / data_csv["p'___(kPa)"]
    data_x = data_csv["p'___(kPa)"]
    data_y = data_csv['q____(kPa)']
    data_z = data_csv['ea_base(%)']
    data_z_riv = data_y / data_x ** (1.6 * 1 /(60 / data_x) ** 0.3)

    fig = go.Figure(data=[go.Scatter3d(x=data_x, y=data_y, z=data_z, marker=dict(
                size=2,    # set color to an array/list of desired values
                opacity=0.5
                )
            ), go.Scatter3d(x=data_x, y=data_y, z=data_z_riv, marker=dict(
                size=2,    # set color to an array/list of desired values
                opacity=0.5,
                color="aqua"
                )
            )
    ])

    fig.update_layout(title='Vs Surface', autosize=False,
                        width=900, height=900,
                        margin=dict(l=65, r=50, b=65, t=90),
                        scene = dict(
                            xaxis = dict(range = [0, 100]),
                            yaxis = dict(range = [-15, 15]),
                            zaxis = dict(range = [-5, 5])
                        )
                    )

    fig.show()

if __name__ == "__main__":
    main()