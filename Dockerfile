from ubuntu:20.04
# install dependencies for the project

RUN apt-get update && apt-get install -y \
    python3 \
    python3-pip

RUN pip3 install --upgrade pip
RUN pip3 install pysam

# copy the project files into the container
COPY . /app
WORKDIR /app

# run the project
CMD ["python3", "main.py"]