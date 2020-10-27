FROM nanome/plugin-env

ENV ARGS=''

ARG CACHEBUST
RUN pip install nanome

COPY . /app
WORKDIR /app

CMD python run.py ${ARGS}
