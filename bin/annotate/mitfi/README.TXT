Version 0.1 of MiTFi
------------------------

The mitochondrial tRNA finder (MiTFi) script is a java script 
using the Infernal package* and several preprocessing steps to find 
mitochondrial tRNA genes. All required covariance models and a file 
defining the NCBI genetic codes** are attached. 


Run the script
------------------------

- You need to have the Infernal package* installed on your system. 
  Please, set the cmsearch call in mitfi_config.txt (default: "cmsearch")

- If you want to run the script on multiple cores (recommended!) you 
  need to define the mpirun/orterun call in mitfi_config.txt 
  (default: "orterun")

- Run the MiTFi script:

   java -jar mitfi.jar mitogenome.fasta 

- Use -h to display all options
   java -jar mitfi.jar -h


We need your feedback! 
------------------------

Please, send your comments, suggestions or questions to
frank@bioinf.uni-leipzig.de



*http://infernal.janelia.org/
**ftp://ftp.ncbi.nih.gov/entrez/misc/data/gc.prt
(Explanation: http://www.ncbi.nlm.nih.gov/Taxonomy/Utils/wprintgc.cgi)


