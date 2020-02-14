FROM continuumio/miniconda3

ENV ARGS=''

COPY . /app
WORKDIR /app

RUN pip install requests
RUN pip install cairosvg
RUN pip install nanome
RUN conda install -c rdkit rdkit

CMD python run.py -a ${ARGS}
