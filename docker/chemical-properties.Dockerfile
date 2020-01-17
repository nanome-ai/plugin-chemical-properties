FROM continuumio/miniconda3

ENV PLUGIN_SERVER=plugins.nanome.ai

COPY . /app
WORKDIR /app

RUN pip install requests
RUN pip install cairosvg
RUN pip install nanome
RUN conda install -c rdkit rdkit

CMD python -m nanome_chemical_properties.ChemicalProperties -a ${PLUGIN_SERVER}