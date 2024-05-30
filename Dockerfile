from ubuntu:20.04
# install dependencies for the project
ARG DEBIAN_FRONTEND=noninteractive
ENV TZ=Etc/UTC

RUN apt-get update
#https://stackoverflow.com/a/73243594

RUN set -ex && \
    apt install -y \
        software-properties-common && \
    add-apt-repository -y ppa:deadsnakes/ppa && \
    apt install -y \
        python3.9 \
        python3.9-distutils \
        python3.9-venv && \
    python3.9 --version && \
    python3.9 -m ensurepip && \
    pip3.9 --version

RUN pip3.9 install --upgrade pip
#RUN pip3 install pysam
# install the project dependencies in requirements.txt
COPY requirements.txt /app/requirements.txt
WORKDIR /app
RUN pip3 install -r requirements.txt
# copy the project files into the container
COPY . /app
WORKDIR /app

# run the project
CMD ["python3.9", "main.py"]