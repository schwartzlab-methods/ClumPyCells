from metadata import *
import plotly.graph_objects as go
import numpy as np
import pandas as pd
import plotly.express as px


def plotInteractiveImage(
    dataFile, colInfo, imageName, area_colName, selected_cellTypes=None, saveName=None
):
    # Unpack column names from colInfo
    imageNumber_colName, x_colName, y_colName, cellType_colName, x_range, y_range = (
        colInfo
    )

    # Filter data for the selected image
    image_data = dataFile[dataFile[imageNumber_colName] == imageName].copy()

    if selected_cellTypes:
        image_data = image_data[image_data[cellType_colName].isin(selected_cellTypes)]

    image_data["radius"] = np.sqrt(image_data[area_colName] / np.pi) * 0.8

    fig = go.Figure()

    unique_cell_types = image_data[cellType_colName].unique()
    colors = px.colors.qualitative.Plotly
    color_map = {
        cell_type: colors[i % len(colors)]
        for i, cell_type in enumerate(unique_cell_types)
    }

    # Adding circles as shapes to the plot
    shapes = []
    for _, row in image_data.iterrows():
        shapes.append(
            {
                "type": "circle",
                "xref": "x",
                "yref": "y",
                "x0": row[x_colName] - row["radius"],
                "y0": row[y_colName] - row["radius"],
                "x1": row[x_colName] + row["radius"],
                "y1": row[y_colName] + row["radius"],
                "line_color": color_map[row[cellType_colName]],
                "fillcolor": color_map[row[cellType_colName]],
                "opacity": 0.5,
            }
        )

    fig.update_layout(
        title="Interactive Plot of Cells",
        xaxis=dict(title="X Coordinate", range=[0, x_range]),
        yaxis=dict(
            title="Y Coordinate", range=[0, y_range], scaleanchor="x", scaleratio=1
        ),
        shapes=shapes,
        showlegend=False,  # Set to True if you want to manually handle legends
        autosize=False,
        width=x_range,
        height=y_range,
    )

    # Save or show the plot
    if saveName:
        fig.write_html(saveName)
    else:
        fig.show()


def plotInteractiveImage3(
    dataFile, colInfo, imageName, area_colName, selected_cellTypes=None, saveName=None
):
    # Unpack column names from colInfo
    imageNumber_colName, x_colName, y_colName, cellType_colName, x_range, y_range = (
        colInfo
    )

    # Filter data for the selected image
    image_data = dataFile[dataFile[imageNumber_colName] == imageName]
    if image_data.empty:
        print("No data found for the specified image name.")
        return

    # Apply filtering for selected cell types if provided
    if selected_cellTypes:
        image_data = image_data[image_data[cellType_colName].isin(selected_cellTypes)]
        if image_data.empty:
            print("No data found for the specified cell types.")
            return

    # Calculate radius from area, ensuring it matches coordinate scale
    image_data["radius"] = np.sqrt(image_data[area_colName] / np.pi) * 0.83

    # Define the figure
    fig = go.Figure()

    # Create a color map for cell types
    unique_cell_types = image_data[cellType_colName].unique()
    colors = px.colors.qualitative.Plotly
    color_map = {
        cell_type: colors[i % len(colors)]
        for i, cell_type in enumerate(unique_cell_types)
    }

    # Adding scatter points for hover information
    for cell_type in unique_cell_types:
        type_data = image_data[image_data[cellType_colName] == cell_type]
        fig.add_trace(
            go.Scatter(
                x=type_data[x_colName],
                y=type_data[y_colName],
                text=[
                    f"Coordinates: ({x}, {y})"
                    for x, y in zip(type_data[x_colName], type_data[y_colName])
                ],
                mode="markers",
                marker=dict(
                    size=type_data["radius"],
                    color=color_map[cell_type],
                    opacity=0.5,
                ),
                name=cell_type,
                hoverinfo="text+name",
            )
        )

    fig.update_layout(
        title="Interactive Plot of Cells",
        xaxis=dict(
            title="X Coordinate", range=[0, x_range], showgrid=False, fixedrange=True
        ),
        yaxis=dict(
            title="Y Coordinate",
            range=[0, y_range],
            showgrid=False,
            fixedrange=True,
            scaleanchor="x",
            scaleratio=1,
        ),
        showlegend=True,
        autosize=False,
        width=800,
        height=800,
    )

    # Save or show the plot
    if saveName:
        fig.write_html(saveName)
    else:
        fig.show()


def plotInteractiveImageDot(
    dataFile: pd.DataFrame,
    dataFile_colNames: list,
    imageName,
    area_colName=None,
    selected_cellTypes: list = None,
    saveName=None,
):
    # Extract column names
    imageName_colName, x_colName, y_colName, cellType_colName, x_range, y_range = (
        dataFile_colNames
    )
    image = dataFile[dataFile[imageName_colName] == imageName]

    if selected_cellTypes:
        image = image[image[cellType_colName].isin(selected_cellTypes)]

    if area_colName:
        image["size"] = np.sqrt(image[area_colName] / np.pi)
    else:
        image["size"] = 5  # Default size if area column not provided

    # Define color scheme with Plotly's qualitative colors
    colors = px.colors.qualitative.Plotly

    # Create the plot
    fig = px.scatter(
        image,
        x=x_colName,
        y=y_colName,
        color=cellType_colName,
        size="size",
        color_discrete_sequence=colors,
        labels={cellType_colName: "Cell Type"},
        title=imageName,
    )

    fig.update_xaxes(range=[-5, x_range])
    fig.update_yaxes(range=[-5, y_range])

    fig.update_layout(
        legend=dict(
            title="Cell Types",
            itemsizing="constant",
            orientation="h",
            yanchor="bottom",
            y=1.02,
            xanchor="right",
            x=1,
        ),
        margin=dict(l=20, r=20, t=40, b=20),
    )

    # Save or show the plot
    if saveName:
        fig.write_html(saveName)
    else:
        fig.show()


aml_metadata = AML_metadata()
imageNumber = 1
if imageNumber > 36:
    dataFile = aml_metadata.nbm_file
else:
    dataFile = aml_metadata.aml_file
dataFile = pd.read_csv(dataFile)
plotInteractiveImage3(
    dataFile,
    aml_metadata.colInfo,
    imageName=imageNumber,
    area_colName="Area",
    selected_cellTypes=None,
    saveName=HOMEDIR + "/Result/images/interactive2.html",
)
