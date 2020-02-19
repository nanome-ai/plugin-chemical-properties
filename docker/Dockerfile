FROM continuumio/miniconda3

ENV ARGS=''

RUN pip install nanome cairosvg
RUN conda install -c rdkit rdkit

COPY . /app
WORKDIR /app

CMD python run.py ${ARGS}
