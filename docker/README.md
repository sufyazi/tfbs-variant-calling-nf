# Building minimal Docker container image for BCFTools

To make the Nextflow pipeline portable and can be run on computing clusters, we need to create a minimal Docker container image to run BCFTools.

### 1. Create a Dockerfile
Create a Dockerfile that would build a minimal Docker image that contains BCFTools and its dependencies. The Dockerfile is located at `docker/` in this repository.

### 2. Build the Docker image
To build the Docker image, run the following command in the terminal, where the Dockerfile is located:

```bash
docker build -t bcftools-minimal .
```

Run `docker images` to check if the Docker image has been built successfully. You should see `bcftools-minimal` in the list of Docker images.

### 3. Test the newly built Docker image
To ensure that the tools installed in the image are working, we can spin up an interactive container of the Docker image and run the tools. Run the following command in the terminal:

```bash
docker run --platform linux/amd64 -p 8080:80 -it bcftools-minimal
```
Once inside the container, just run `bcftools` to check if the tool is working. You should see the help message of `bcftools` printed in the terminal.

<small>**Note**: the `--platform linux/amd64` flag is used to specify the platform of the Docker image, because we built the container image for Linux machines. This flag is only recommended when building and running Docker images on Apple Silicon (Mx) Macs. If you are using an Intel-based Mac, you can omit this flag.</small>

### 4. Save the Docker image to a file
If everything works as expected, we can save the Docker image to a file. 
    
Run the following command in the terminal:

```bash
docker save -o bcftools-minimal.tar bcftools-minimal
```

This will save the Docker image to a file named `bcftools-minimal.tar` in the current directory.

### 4a. (Optional) Push the Docker image to Docker Hub
If you want to share the Docker image with others, you can push the Docker image to Docker Hub. First, you need to tag the Docker image with your Docker Hub username and the name of the repository. Run the following command in the terminal:

```bash
docker tag bcftools-minimal:tagname DOCKER-USERNAME/bcftools-minimal:tagname
```

Make sure that 1) you already have a Docker Hub account, 2) you have logged in to Docker Hub using `docker login` before running the command above, and 3) you have already created an empty repository on Docker Hub to tag and push the local image to. Replace `DOCKER-USERNAME` accordingly and `tagname` with your chosen tag (or omit `tagname` to use `latest` as the default tag).

Now push the Docker image to Docker Hub:

```bash
docker push DOCKER-USERNAME/bcftools-minimal:tagname
```
<small>**Note**: replace `tagname` with your chosen tag. This helps organize your Docker Hub repository. If you don't specify a tag, Docker will use `latest` as the default tag.</small>


### 5. Convert the Docker image to a Singularity image
As most computing clusters do not support Docker due to [security concerns](https://arctraining.github.io/hpc2-software/course/containers.html), we need to convert the Docker image to a Singularity image. To do this, we need to use Singularity. As we are already using containerization tools, there is no need to bother installing Singularity locally. Instead, we can use the Singularity image from the [BioHPC repository](https://hub.docker.com/r/biohpc/sindocker) and pull a Singularity image from there.

Once done, we can use Docker to run an instance of the Singularity image, and then from within this container, we can use the `singularity` command-line tool to convert the Docker image we just created into a Singularity image (`.sif`). 

With `docker run` we can automate all the steps here, by mounting the current directory where our image lives, then pulling the BioHPC Singularity image, running it, and then running the conversion inside the container. Run the following command in the terminal:

```bash
docker run --platform linux/amd64 --rm -v $(pwd):/data biohpc/sindocker singularity build --force bcftools-minimal.sif docker-archive://bcftools-minimal.tar 
```

### 6.  Test the Singularity image
Next, transfer the `.sif` file to your favourite HPC cluster platform where Singularity has presumably already been installed as a loadable module. Load up Singularity on the cluster, and change the `.sif` file permission before testing the image.

```bash
chmod a+x bcftools-minimal.sif
```

Spin up an interactive container of the Singularity image and run the tool. 

```bash
singularity shell bcftools-minimal.sif
```

Alternatively, run the image directly without entering the shell:

```bash
./bcftools-minimal.sif bcftools --help
```

If everything is done correctly, you should see the help message of `bcftools` printed in the terminal. If the command returned a `command not found` error, it is most likely that there is an issue with the setting of the `$PATH` in the image. This should not happen, but if it did, consider exporting the `$PATH` inside the running container to include the directory where `bcftools` is located.

```bash
export PATH="/opt/conda/bin:/opt/conda/condabin":$PATH
```

### 6a. (Optional) Pull the Docker Image directly from Docker Hub using Singularity
If you have pushed the Docker image to Docker Hub, you can pull the Docker image directly from Docker Hub using Singularity, which will pull the image and convert it into SIF format in one step. Run the following command in the terminal:

```bash
singularity shell docker://DOCKER-USERNAME/bcftools-minimal:tagname
```

A container terminal should start and running `bcftools` should return the help message of `bcftools` printed in the terminal.

___

### Final Notes
The Docker image is now ready to be used in the Nextflow pipeline. Both Singularity and Docker images are supported by Nextflow, so you can choose to use either one of them.


