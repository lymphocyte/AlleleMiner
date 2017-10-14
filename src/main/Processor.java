package main;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.sql.SQLException;
import java.util.ArrayList;
import java.util.Properties;

import algorithm.AlleleDetection;
import input.ReadVCFfile;
import tools.QualityController;
import tools.genmoeSampleFileParser;
import input.ReadGene;

public class Processor {

    private String final_result = "";
    Properties props = new Properties();
    private boolean testDatabase;
    private int total_haplotypes = 0;

    public void Process(String vcf_file_path, String gene_position_file_path, String outputpath, String vcfmode,
	    String threshold) throws FileNotFoundException, IOException, SQLException {

	props.load(Run.class.getResourceAsStream("/src/resource/properties.properties")); //when package
	//props.load(Run.class.getResourceAsStream("/resource/properties.properties")); //when test
	this.props.setProperty("vcf_path", vcf_file_path);
	this.props.setProperty("vcf_source", vcfmode);
	this.props.setProperty("allele_threshold", threshold);
	this.testDatabase = Boolean.valueOf(this.props.getProperty("testDatabase"));
	FileWriter fw = new FileWriter(new File(outputpath));
	/*
	 * Step 1. Get the gene sequence of the interesting gene. User type in
	 * the gene they are interested in. There are two types of typing from
	 * user:1. gene name. Then we search from gene database to find the
	 * chormosome and location. 2. chormosome and location. then we directly
	 * find the gene in the database.
	 */
	// ReadGene reader = new ReadGene("C2","genename",props);
	// ReadGene reader = new
	// ReadGene("chr7,141998851,141999700","location",props);
	// ReadGene reader = new
	// ReadGene("chr7,141998851,142510972","location",props);
	// ReadGene reader = new
	// ReadGene("chr14,106032614,107288051","location",props);
	// ReadGene reader = new
	// ReadGene("chr14,106494134,106494577","location",props);
	// 106494135 106494597 106405611|106406108
	// BufferedReader br = new BufferedReader(new
	// FileReader("IGH_postions"));
	BufferedReader br = new BufferedReader(new FileReader(gene_position_file_path));
	String postions = br.readLine();

	while (postions != null) {
	    String chr = postions.split(",")[0];
	    if (chr.contains("HG")) {
		postions = br.readLine();
		continue;
	    }
	    /*
	     * For every gene, we need to create a folder, in every folder, it
	     * has the alleles from each haplotype.
	     */
	    String gene_name = postions.split(",")[1];
	    String folder_path = "/Original_haplotypes/";
	    File folder = new File(folder_path);
	    folder.mkdir();
	    System.out.println(gene_name);
	    String start_index = postions.split(",")[2];
	    String end_index = postions.split(",")[3];
	    // String start_index = "142270879";
	    // String end_index = "142271415";
	    ReadGene reader = new ReadGene(chr + "," + start_index + "," + end_index + "", "location", props,
		    postions.split(",")[4]);
	    /*
	     * Step 2. Read The VCF file, which only contains the SNPs of the
	     * target gene.
	     */
	    System.out.println(reader.GetDNASeq());
	    ReadVCFfile vcf_file = new ReadVCFfile(reader.GetStartIndex(), reader.GetEndIndex(), reader.GetChr(),
		    props);
	    ArrayList<String> indexed_snp_results = vcf_file.getIndexedResults();
	    String[] sample_IDs = vcf_file.getSampleIDs();
	    total_haplotypes = sample_IDs.length*2;
	    genmoeSampleFileParser gepar = new genmoeSampleFileParser(
		    "/Users/yuyaxuan/Documents/workspace/AlleleMiner/src/resource/sample_info.csv");
	    if (this.props.getProperty("testDatabase").contains("true")) {
		gepar.read_sample_file();
	    }
	    /*
	     * Step 3.
	     */
	    AlleleDetection ald = new AlleleDetection(reader.GetDNASeq(), indexed_snp_results, reader.GetStartIndex(),
		    reader.GetEndIndex(), gene_name, "", sample_IDs, gepar.getSample_map(), this.testDatabase,
		    folder_path, reader.get_reverse_info());
	    this.final_result = ald.getFinal_result();
	    fw.write(this.final_result);
	    postions = br.readLine();
	}

	// System.out.println(this.final_result);

	fw.close();
    }

    public int GetTotalHaplotypes() {
	return total_haplotypes;
    }
}
