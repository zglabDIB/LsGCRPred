<h1>LsGCRPred(LSNP-Gene Interaction Based Cancer Risk Prediction Model)</h1>

<h2>Introduction</h2>

<p align="center">
  <img src="assets/arch_git.png" alt="Project Screenshot" width="600" height="700">
</p>

<b>Breast and Ovarian cancer risk prediction models</b> consolidated in a single tool. It looks at lncRNA target gene expression differences between cancer and normal tissues, along with the presence or absence of key variants. Using <b>9 lncRNA-associated SNPs(2 specific to Breast Cancer and 7 specific to Ovarian cancer)</b>, it helps distinguish between high-risk and non-cancerous states.

<h2>Pre-requisite</h2>
The input provided by user should be paired-end RNA-Seq data corresponding to Breast or Ovarian tissue of epithelial origin.<br/>
The BAM pre-processing tool <b>Opossum</b> is required for preparing RNA-seq BAM files prior to variant calling with <b>Platypus</b>. Both Opossum and Platypus are based on Python 2.7 and is not compatible with Python 3.<br/> 
A dedicated Conda environment for Python 2.7 with several dependencies is therefore required to run <b>Workflow 1</b>. Please follow the environment setup instructions provided below.<br/>
<b>##Note:</b> GATK pipeline can also be utilized. The results highly overlap while Platypus and Opossum are computationally much faster.

<h2>Input Data format</h2>
<b>RNA-seq Data:</b> Paired-end RNA-seq data must be pre-processed to identify SNPs and generate FPKM values, which will be used as input for the prediction model. The following steps outline the required pre-processing workflow.

<h3>A. Workflow 1: Generation of Gene expression profile and Variant file</h3>
<ol type="1">
    <li>Python programming Environment setup: Download <b>Anaconda</b> or <b>Miniconda</b> distribution from their official website</li>
    <li>
        Create and activate the conda environment and install dependencies
        <pre>conda create -n LsGCRPred_VARIANT_EXP -c conda-forge -c bioconda python=2.7 cython=0.29 pysam=0.15 samtools=1.9 htslib hisat2 star stringtie fastp bcftools vcftools</pre>
        #Activate the conda environment using
        <pre>conda activate LsGCRPred_VARIANT_EXP</pre>
        #OR
        <pre>source activate LsGCRPred_VARIANT_EXP</pre>
    </li>
    <li>
        Move to the <b>Workflow1</b> folder within the downloaded <b>LsGCRPred</b> directory using
        <pre>cd path_to_folder/LsGCRPred/Workflow1</pre>
    </li>
    <li>
        Download <b>Opossum</b> from github
        <pre>git clone https://github.com/CGGOxford/Opossum.git</pre>
    </li>
    <li>
        Download <b>Platypus</b> Variant caller from github and build it.<br/>
        #Platypus requires htslib and Cython (installed in Step 2). Once installation is complete, Platypus.py will be located in ./Platypus/bin/. For any installation issues, refer to the official Platypus GitHub documentation, as setup may vary across systems.
        <pre>
            git clone https://github.com/andyrimmer/Platypus.git
            export CFLAGS="-I$CONDA_PREFIX/include"
            export LDFLAGS="-L$CONDA_PREFIX/lib"
            export CPATH=$CONDA_PREFIX/include
            cd Platypus
            make
        </pre>
        #Leave the Platypus folder and return to the <b>Workflow1</b> folder
        <pre>cd ../</pre>
    </li>    
    <li>
        Enter the folder <b>reference</b> within <b>Workflow1</b>
        <pre>cd reference</pre>
    </li>
    <li>
        Download <b>hg38genome fasta and v49 gtf annotation</b> filesfrom https://www.gencodegenes.org/human/ in this folder<br/>
        #fasta
        <pre>
            wget https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_49/GRCh38.primary_assembly.genome.fa.gz
            gunzipGRCh38.primary_assembly.genome.fa.gz
        </pre>
        #gtf<br/>
        <b>#Note:</b>The current models were trained using transcript IDs from GENCODE v49. If future GTF releases modify or remove any of these transcripts, the model may not function correctly for the affected SNPs. We therefore recommend using GENCODE v49 from the archive. For further assistance, please contact the authors.
        <pre>
            wget https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_49/gencode.v49.annotation.gtf.gz
            gunzip gencode.v49.annotation.gtf.gz
        </pre>
    </li>
    <li>
        Generate genome fasta index files<br/>
        <b>#Note:</b>We use 8 threads for this workflow by default; adjust this based on your available resources. Users can choose either HISAT2 or STAR as the aligner.
        <ol type="a">
            <li>
                If using Hisat2
                <pre>hisat2-build -p 8 GRCh38.primary_assembly.genome.faGRCh38.primary_assembly.genome</pre>
            </li>
            <li>
                If using STAR
                <pre>STAR --runThreadN 8 --runMode genomeGenerate --genomeSAindexNbases 14 --genomeDir ./STAR --genomeFastaFiles GRCh38.primary_assembly.genome.fa --sjdbGTFfile gencode.v49.annotation.gtf</pre>
            </li>
        </ol>    
    </li>    
    <li>
        Download the compressed dbsnp <b>reference variant annotation and tab indexed files</b> from ddbj
        <pre>
            wget https://ddbj.nig.ac.jp/public/public-human-genomes/GRCh38/fasta/dbsnp_146.hg38.vcf.gz
            wget https://ddbj.nig.ac.jp/public/public-human-genomes/GRCh38/fasta/dbsnp_146.hg38.vcf.gz.tbi
        </pre>
        #Leave the folder
        <pre>cd ../</pre>
    </li>
</ol>

Perform the below steps to generate <b>TEST_FPKM.txt</b> (transcript expression file) and <b>TEST_{breast/ovary}_selected_SNP.txt</b>(variant file)
<b>#Note:</b>The word “TEST” will be replaced with your sample name.<br/>

<ol start="10">
    <li>
        Enter the folder <b>samples</b> within <b>Workflow1</b>
        <pre>cd samples</pre>
    </li>
    <li>
        Place paired-end Breast/Ovary FASTQ files in this directory. Supported file formats: .fastq, .fastq.gz, .fq, .fq.gz.<br/>
        Example: <b>TEST_1.fastq.gz, TEST_2.fastq.gz.</b> Multiple samples can be included.
    </li>
    <li>
        Leave the directory and return to the <b>Workflow1</b> folder
        <pre>cd ../</pre>
    </li>
    <li>
        Edit file <b>sample_ids.txt</b> within <b>Workflow1</b> folder and write one sample name per line (ex. Replace <b>TEST</b> with your sample name)
    </li>
    <li>
        Run the bash script from Workflow1 folder<br/>
        <ol type="a">
            <li>
                If using Hisat2<br/>
                #for breast
                <pre>bash exp_var.sh -g ./reference -s ./sample_ids.txt -a hisat2 -i ./samples -o output_variant -t 8 -c breast</pre>
                #for ovary
                <pre>bash exp_var.sh -g ./reference -s ./sample_ids.txt -a hisat2 -i ./samples -o output_variant -t 8 -c ovary</pre>
            </li>
            <li>
                If using Star<br/>
                #for breast
                <pre>bash exp_var.sh -g ./reference -s ./sample_ids.txt -a star -i ./samples -o output_variant -t 8 -c breast</pre>
                #for ovary
                <pre>bash exp_var.sh -g ./reference -s ./sample_ids.txt -a star -i ./samples -o output_variant -t 8 -c ovary</pre>
                <b>#Note:</b> default cpu usage is 8 unless mentioned otherwise. For detailed usage, check <b>bash exp_var.sh -h</b>
                #Make sure you are in <b>Workflow1</b> folder where Opossum and Platypus both are installed.
            </li>
        </ol>    
    </li>    
    <li>
        Go to the output folder
        <pre>cd output_variant/TEST</pre>
    </li>
    <li>
        Check for the presence of
        <ol type="a">
            <li>
                Transcript expression file TEST_FPKM.txt
            </li>
            <li>
                Selected variants from file <b>TEST_{breast/ovary}_selected_SNP.txt</b>     
            </li>
        </ol>    
    </li>
</ol>

<ol start="17">
    <li>
        Deactivate this conda environment
        <pre>Conda deactivate</pre>
    </li>
</ol>    


<h3>SelectedBreast SNP:</h3>
rs2366152<br/>
rs7091441<br/>
<h3>SelectedOvarian SNP:</h3>
rs932501<br/>
rs2072588<br/>
rs7318592<br/>
rs9506960<br/>
rs9510420<br/>
rs12583808<br/>
rs60135126<br/><br/>
<b>#Note:</b> If no SNP is detected, the file would be empty, in that case, no need to proceed to <b>Workflow2</b>.


<h3>B.	Workflow 2: Executing the prediction model</h3>

Create a conda environment named <b>LsGCRPred_model</b> using one of the commands below. Both Python 3.8 and Python 3.9 options are provided, use whichever is compatible with your server environment.

<pre>
    conda create -n LsGCRPred_model -c conda-forge python=3.9 numpy=1.21 pandas=1.3 scikit-learn=1.1.3 pickleshare=0.7.5 
</pre>

#Activate the conda environment using 
<pre>conda activate LsGCRPred_model</pre>
#OR
<pre>source activate LsGCRPred_model</pre>

#OR
<b>For installing the packages involving python 3.8 following procedure is required to be followed:</b>

<ul>
    <li>
        Create a conda environment named viz. LsGCRPred_model using the following command in order to install python 3.8
        <pre>conda create -n LsGCRPred_model -c bioconda -c conda-forge python=3.8 -y</pre>
    </li>
    <li>
        Activate the conda environment using
        <pre>conda activate LsGCRPred_model</pre>
        #OR
        <pre>source activate LsGCRPred_model</pre>
    </li>
    <li>
        Download the necessary python packages from the requirement.txt file provided here using
        <pre>pip  install -r requirement.txt</pre>
    </li>
</ul>

<ol start="1">
    <li>
        <b>Enter the main LsGCRPred directory using</b>
        <pre>cd path_to_folder/LsGCRPred/</pre>
    </li>
    <li>
        <b>Program Execution Procedure:</b><br/><br/>
        Two separate Python scripts are provided for model execution: <b>breast_snps_calling.py</b> and <b>ovary_snps_calling.py</b>, corresponding to breast and ovarian tissue, respectively.<br/>
        The input expression file (generated in Step 14 of Workflow 1 as <b>TEST_FPKM.txt</b>) should be placed in the <b>input_files</b> directory within both the breast and ovary folders.<br/>
        Rename the expression file as follows:<br/>
        <b>input_expression_breast.txt</b> (for breast)<br/>
        <b>input_expression_ovary.txt</b> (for ovary)<br/>
        The file must contain two columns with the headers: <b>t_name</b> and <b>FPKM</b>.
        Example input files for both breast and ovary are provided in the example folder.  
    </li>
    <ol type="a">
        <li>
            To execute the Breast tissue sample, enter the following:
            <pre>python breast_snps_calling.py</pre>
        </li>
        <li>
            To execute theOvarian tissue sample, enter the following:
            <pre>python ovary_snps_calling.py</pre>
        </li>
    </ol>
    During execution, you will be prompted to specify whether SNPs are present <b>(Y or N)</b>. Based on your input, the corresponding output files will be generated in the output_files directory for the selected tissue.        
</ol>    

<h2>Results</h2>
<ul>
    <li>
        The output file will be named <b>snp_id_pred_output.txt</b> (e.g., <b>rs2366152_pred_output.txt</b> for a breast sample).
    </li>
    <li>
        The output includes a <b>pred_proba</b> column, which represents the predicted probability of cancer risk.
    </li>
    <li>
        If <b>(pred_proba*100) >=50</b>,the prediction is considered positive (i.e., a probable risk of cancer).
    </li>
</ul>
