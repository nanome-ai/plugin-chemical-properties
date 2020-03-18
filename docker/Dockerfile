FROM continuumio/miniconda3

ENV ARGS=''

RUN conda install -c rdkit rdkit
RUN pip install cairosvg

ARG CACHEBUST
RUN pip install nanome

COPY . /app
WORKDIR /app

CMD python run.py ${ARGS}
