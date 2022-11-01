FROM continuumio/miniconda3

RUN apt-get update && apt-get install --quiet --yes --no-install-recommends libgl1 && apt-get clean && rm -rf /var/lib/apt/lists/* && mkdir /install
COPY environment.yml setup.py vsflow README.md /install/
COPY vslib /install/vslib
RUN conda env create --quiet --file /install/environment.yml && conda clean --all --yes --quiet
RUN bash -c "source /opt/conda/bin/activate vsflow && pip install /install" && rm -rf /install && sed -ie 's/base/vsflow/' /root/.bashrc
