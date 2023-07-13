# MitoZ 3
<img alt="GitHub release (latest SemVer)" src="https://img.shields.io/github/v/release/linzhi2013/mitoz?label=Latest%20release"> <img alt="GitHub all releases" src="https://img.shields.io/github/downloads/linzhi2013/MitoZ/total?label=Github%20downloads">  <img alt="GitHub" src="https://img.shields.io/github/license/linzhi2013/mitoz?label=License">    
[![docker build](https://img.shields.io/badge/docker%20build-passing-brightgreen)](https://hub.docker.com/r/guanliangmeng/mitoz/tags) [![docker latest version](https://img.shields.io/docker/v/guanliangmeng/mitoz)](https://hub.docker.com/r/guanliangmeng/mitoz/tags) [![docker pulls](https://img.shields.io/docker/pulls/guanliangmeng/mitoz?style=flat)](https://hub.docker.com/r/guanliangmeng/mitoz/tags) 
[![singularity image](https://img.shields.io/badge/Singularity%20build-passing-brightgreen)](https://github.com/linzhi2013/MitoZ/wiki/Installation#3-apptainersingularity)      
[![Conda-pack](https://img.shields.io/badge/conda--pack-passing-brightgreen)](https://github.com/linzhi2013/MitoZ/wiki/Installation#4-conda-pack)  [![install with bioconda](https://img.shields.io/badge/install%20with-bioconda-brightgreen.svg?style=flat)](https://github.com/linzhi2013/MitoZ/wiki/Installation#5-conda)   [<a href="https://anaconda.org/bioconda/mitoz"> <img src="https://anaconda.org/bioconda/mitoz/badges/version.svg" /> </a>](https://anaconda.org/bioconda/mitoz/badges/version.svg)   [![Anaconda-Server Badge](https://img.shields.io/conda/dn/bioconda/mitoz)](https://anaconda.org/bioconda/mitoz)   


THE SOFTWARE IS PROVIDED "AS IS" AND THE AUTHOR DISCLAIMS ALL WARRANTIES WITH
REGARD TO THIS SOFTWARE INCLUDING ALL IMPLIED WARRANTIES OF MERCHANTABILITY
AND FITNESS. IN NO EVENT SHALL THE AUTHOR BE LIABLE FOR ANY SPECIAL, DIRECT,
INDIRECT, OR CONSEQUENTIAL DAMAGES OR ANY DAMAGES WHATSOEVER RESULTING FROM
LOSS OF USE, DATA OR PROFITS, WHETHER IN AN ACTION OF CONTRACT, NEGLIGENCE
OR OTHER TORTIOUS ACTION, ARISING OUT OF OR IN CONNECTION WITH THE USE OR
PERFORMANCE OF THIS SOFTWARE.

**WHEN YOU ADAPT (PART OF) THE SOFTWARE FOR YOUR USE CASES, THE AUTHOR AND
THE SOFTWARE MUST BE EXPLICITLY CREDITED IN YOUR PUBLICATIONS AND SOFTWARE,
AND YOU SHOULD ASK THE USERS OF YOUR SOFTWARE TO CITE THE SOFTWARE IN
THEIR PUBLICATIONS. IN A WORD, 请讲武德.**

**About:**
- **MitoZ provides a "one-click" solution to get annotated mitogenomes from raw data fastq files.**

**News**
- **(April-20-2023) [Docker](https://github.com/linzhi2013/MitoZ/wiki/Installation#1-docker) and [Singularity](https://github.com/linzhi2013/MitoZ/wiki/Installation#3-apptainersingularity) versions of MitoZ 3.6 come out now. Both were tested on Ubuntu 20.04.4 LTS.**

- **(April-19-2023) On the installation problem (https://github.com/linzhi2013/MitoZ/issues/188), hopefully now it is fixed. Please let me know if it is not.**

- **(April-14-2023) MitoZ 3.6 is just released (https://github.com/linzhi2013/MitoZ/releases/tag/3.6), fixed some bugs in MitoZ 3.5. It is recommended to upgrade to this version!** You can install it via **[conda-pack](https://github.com/linzhi2013/MitoZ/wiki/Installation#4-conda-pack)** (firstly recommended if the conda way does not work for you), [conda](https://github.com/linzhi2013/MitoZ/wiki/Installation#5-conda) and [source code](https://github.com/linzhi2013/MitoZ/wiki/Installation#6-source-codes).

**See**
- **Installation: https://github.com/linzhi2013/MitoZ/wiki/Installation**.
- **MAKE SURE that you do a [test run](https://github.com/linzhi2013/MitoZ/wiki/Installation#9-running-the-test-dataset) using provided test dataset before running your own samples!**. 
- **Tutorial: https://github.com/linzhi2013/MitoZ/wiki/Tutorial (Recommended if you are new to MitoZ!)**
- **Documentation: https://github.com/linzhi2013/MitoZ/wiki** and the [HMTL version](https://github.com/linzhi2013/MitoZ/releases/download/3.6/MitoZ_manual_v3.6.html). The HTML version may not be update-to-date.
- **Latest release: https://github.com/linzhi2013/MitoZ/releases/**


**Bugs and Questions**
- Have a look at https://github.com/linzhi2013/MitoZ/issues and https://github.com/linzhi2013/MitoZ/wiki/Known-issues for known bugs or issues.

- **Please try the latest version first if you find some bugs in the old versions**
	- to do that, you should specify the version of MitoZ when you use the `mamba/conda` command (please refer to the [installation instruction](https://github.com/linzhi2013/MitoZ/wiki/Installation)), as I found out that many people still download the older versions.

- **Known bugs for MitoZ 3.5 (April-13-2023): (1) If your default shell is not bash, you can run into the missing annotation of tRNA genes (see https://github.com/linzhi2013/MitoZ/issues/187). Please change the default shell to bash before using MitoZ 3.5! (2) In MitoZ 3.5, I mistakenly used a `cmsearch` binary for Mac OS for Linux platform, which leads to the problem of failing to annotate any tRNA genes. Please check https://github.com/linzhi2013/MitoZ/issues/187 for the current solution.** 

- check which shell you are using:
    ```
    $ echo "$SHELL"
    ```


- _I have been updating the documentation (wiki) from time to time, so it may be good for you to check the documentation again every after some time._

- In case there are still bugs in the latest version, please **firstly search https://github.com/linzhi2013/MitoZ/issues** and the Wiki to check whether similar questions have been raised by other users. If no related issues and answers are found, then please raise a new issue. Thank you!

- Any feedbacks are wellcome!


# Citations
- Meng G, Li Y, Yang C, Liu S. MitoZ: a toolkit for animal mitochondrial genome assembly, annotation and visualization. Nucleic acids research. 2019 Jun 20;47(11):e63-. https://doi.org/10.1093/nar/gkz173
- Additionally, please cite the related software invoked by MitoZ: https://github.com/linzhi2013/MitoZ/wiki/Citations
