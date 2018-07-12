# Using TOPMed's Docker Image for the TOPMed's workshop
Users have two choices for doing the exercises in the workshop:
1. Logging into an RStudio session on a remote computer.
2. Run a TOPMed docker image on your local computer.

This project describes how to do the exercises using the docker image.
## Prerequisites ##
The following software is required on your local computer or laptop:
- [ ] git is required
- [ ] Python 2.7.11+ is required
- [ ] Docker engine (aka <i>docker</i>)  

  *  For Mac computers, installing <i>docker</i> is fairly straight forward.  Following the Docker instructions at the url https://docs.docker.com/docker-for-mac/install/

  * For Windows 10 computers, installing docker requires enabling <i>Hyper-V</i> in Windows and enabling virtualization in the bios. Follow the Docker instructions at the url: https://docs.docker.com/docker-for-windows/install/
## General Steps ##
The general steps for running RStudio using the docker images are as follows:
1. Download the docker image
2. Create a working directory and <i>cd</i> to it
2. Download the docker helper scripts
3. Running RStudio server in the docker image
4. Connect to the RStudio server using the local browser


#### Download the Docker Image ####
After <i>docker</i> has been installed, download the TOPMed RStudio docker image <kbd>uwgac/topmed-RStudio</kbd>.  For example:
```
docker pull uwgac/topmed-RStudio
```
This step only has to be done once.  The image will be stored persistently on the computer.

#### Create Working Directory ####
Create a working directory and <i>cd</i> to it.  For example,
```
mkdir ~/workshop_2018
cd ~/workshop_2018
```
This step only has to be done once.

#### Download the docker helper functions ####
Using <i>git</i>,  download the docker help functions and create an alias to the helper function <i>Rstudio_docker.py</i>.  For example,
```
mkdir ~/workshop_2018
cd ~/workshop_2018
git clone https://github.com/uw-gac/docker_helpers
alias rs_docker='~/workshop_2018/docker_helpers/Rstudio_docker.py'
```

This step also only has to be done once.

#### Running RStudio ####
Run RStudio server within a docker container in the background (or detached) using the helper function <i>Rstudio_docker.py</i> (e.g., using the previously defined alias).  This enables a browser on a local computer to start an RStudio session.

The RStudio session has the following characteristics and attributes:
1. The user name and password is <i>rstudio</i>
2. The current working directory of the RStudio session is <i>/home/rstudio</i>
3. By default the current working directory in the RStudio session is mapped to the current working directory of the local computer

The general command syntax is:
```{r}
rs_docker [options]
```
Execute <i>rs_docker --help</i> for more details.

<i>Examples</i>
1. <i>Example 1</i>
```{r}
mkdir ~/workshop_2018
cd ~/workshop_2018
git clone https://github.com/uw-gac/docker_helpers
alias rs_docker='~/workshop_2018/docker_helpers/Rstudio_docker.py'
rs_docker
```
The example above runs RStudio server in the background in the container named <i>rstudio</i>. The current working directory, <i>~/workshop_2018/</i>, will be mapped to <i>/home/rstudio</i> in the docker container.

In the local computer's browser, the following URL will create an RStudio session:
```{r}
localhost:8787/
```
The current working directory in the R session is <i>/home/rstudio</i> which is also mapped to the local computer's directory <i>~/workshop_2018/</i>.

2. <i>Example 2</i>
```{r}
rs_docker -i 192.168.1.4 -p 8080
```
The example above runs RStudio server in the background in the container named <i>rstudio</i>.  In a local computer's browser, the following URL will create an RStudio session:
```{r}
192.168.1.4:8080/
```
3. <i>Example 3</i>
```{r}
rs_docker -C kill
```
The example above kills the background docker container named <i>RStudio</i>
