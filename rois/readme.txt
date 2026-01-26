This folder cointains atlases used for defining ROIs and for anatomical
labeling. These files are not needed for rsfc analyses. The atlases for
anatomical labeling need to be located in a subfolder called 'atlas'
which needs to contain the following files (which can be obtained from
various sources):

aal2.nii (Rolls et al., 2015)
aal2.txt
AICHA.nii (Joliot et al., 2015)
AICHA.txt
BNA.nii (Fan et al., 2016)
BNA.txt
brodmann.nii (from MRICron)
brodmann.txt
cerebellum.nii (Diedrichsen et al., 2009)
cerebellum.txt
HarvardOxford-Cortical.nii (from FSL)
HarvardOxford-Cortical.txt
HarvardOxford-Subcortical.nii (from FSL)
HarvardOxford-Subcortical.txt
lpba40.nii (Shattuck et al., 2008)
lpba40.txt
natbrainlab.nii (de Schotten et al., 2011)
natbrainlab.txt
subcortical.nii (Pauli et al., 2017)
subcortical.txt
thalamus.nii (Najdenovska et al., 2018)
thalamus.txt


The txt files need to contain the labels corresponding to the nii files
in the following format:

1 Label1
2 Label2
3 Label3
...
