package algorithm; //file is included in package algorithm

//include java built-in classes
import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashMap;
import java.util.Iterator;
import java.util.Random;

//another package is included
import tools.ReverseCompliment;

public class AlleleDetection { //main class of the file
    private String ref_seq;
    // private String snp_raw_result;
    private ArrayList<String> indexed_snps;
    private ArrayList<Haplotype> Haplotype_factory = new ArrayList<Haplotype>(); //an ArrayList of Haplotype which is another in algorithm package
    private int[] snp_positions;
    private boolean testDatabase; // This will decide if it's in the test mode,
				 // which will randomly reduce 10 percent of the
				 // amount of Haplotypes.  
    private double percentage_of_alleles_to_reduce = 0.1;
    private int sample_size;
    private int info_length = 9;
    private int snp_num = 0;
    private int start_index;
    private int end_index;
    private String genename;
    private String final_result = "";
    private String[] sample_ids;
    private String folder_path;
    private String rev_info = null;
    private HashMap<String, String> sample_map = new HashMap<String, String>(); //hash map with key and value both are strings

    //parameterized constructor
    public AlleleDetection(String ref_gene_seq, ArrayList<String> snp_indexed_result, int start_index, int end_index,
	    String genename, String result, String[] sample_IDs, HashMap<String, String> sample_map,
	    boolean testDatabase, String folder_path, String rev_info) throws IOException { //IO Excetion is for File handling classes
	this.final_result = result; //this is used to refer the data members of the class
	this.folder_path = folder_path;
	this.genename = genename;
	this.ref_seq = ref_gene_seq;
	this.indexed_snps = snp_indexed_result;
	this.sample_ids = sample_IDs;
	this.sample_map = sample_map;
	this.rev_info = rev_info;
	this.testDatabase = testDatabase;
	System.out.println(this.genename + " Haplotypes between " + start_index + " : " + end_index);
	this.final_result += this.genename + " Haplotypes between " + start_index + " : " + end_index + " "
		+ this.rev_info + System.lineSeparator();
	if (snp_indexed_result.size() == 0) {
	    String tmp_seq = "Allele" + 1 + " " + this.ref_seq + " count is 5008 ***Main allele***"
		    + System.lineSeparator(); //line seperator for new line
	    System.out.println(tmp_seq);
	    this.final_result += tmp_seq;
	    return;
	}
	this.snp_positions = new int[this.indexed_snps.size()];
	this.sample_size = this.indexed_snps.get(0).split("\t").length - info_length;
	String[] test = this.indexed_snps.get(0).split("\t");
	this.start_index = start_index;
	this.end_index = end_index;
	this.snp_num = this.indexed_snps.size();
	Process();
    }

    private void Process() throws IOException {
	/*
	 * Create 2*sample_size haplotypes.
	 */

	Haplotype[] snp_haplotypes = new Haplotype[this.sample_size * 2];
	/*
	 * Initialize Haplotypes
	 */

	for (int i = 0; i < snp_haplotypes.length; i++) {
	    Haplotype hap = new Haplotype(this.ref_seq, this.start_index, this.end_index, this.genename);
	    // hap
	    hap.setID(this.sample_ids[i / 2]);
	    if (this.sample_map.size() != 0) {
		hap.setOrigin(this.sample_map.get(hap.getID()));
	    }
	    snp_haplotypes[i] = hap;
	}

	for (int i = 0; i < this.indexed_snps.size(); i++) {
	    String[] temp_content = this.indexed_snps.get(i).split("\t");
	    this.snp_positions[i] = Integer.parseInt(temp_content[1]) - this.start_index;
	    if (this.ref_seq.charAt(this.snp_positions[i]) != temp_content[3].charAt(0)) {
		    //charAt: function used to return the character from specific index
		String ref_nt = this.ref_seq.charAt(this.snp_positions[i]) + "";
		ReverseCompliment RevC = new ReverseCompliment(temp_content[3].charAt(0) + "");
		String vcf_nt = RevC.getReverseCompliment(); //according to my rnd, reverse compliment for DNA is described at the end of file
		if (vcf_nt.charAt(0) == ref_nt.charAt(0)) {
		    // If the VCF file is reversed, then we reverse it back to
		    // the right order.
		    temp_content[3] = vcf_nt;
		    temp_content[4] = new ReverseCompliment(temp_content[4]).getReverseCompliment();

		} else {
		    System.out.println("Something wrong with the index!");
		    System.out.println(temp_content[2] + "-" + temp_content[1] + "-" + temp_content[3] + "==>"
			    + temp_content[4] + "-" + this.ref_seq.charAt(this.snp_positions[i]) + "");
		    System.exit(0); //terminate the program
		    continue; //again go for next iteration
		}

	    }
	    /*
	     * Put 0|1 structures of each SNP in an array.
	     */
	    int[] snp_map = new int[this.sample_size * 2];
	    int tmp_counter = 0;
	    for (int m = this.info_length; m < snp_map.length / 2; m++) { //length: size of snp_map
		if (m >= temp_content.length) {
		    System.out.println("out of border");
		    continue;
		}

		// temp_str is 1|0
		String temp_str = temp_content[m];
		if (temp_str.length() == 3) {
		    snp_map[tmp_counter] = Integer.parseInt(temp_content[m].split("\\|")[0]);//parseInt: convert string to integer
		    tmp_counter++;
		    snp_map[tmp_counter] = Integer.parseInt(temp_content[m].split("\\|")[1]);
		    tmp_counter++;
		} else {

		    snp_map[tmp_counter] = 0;
		    tmp_counter++;

		    snp_map[tmp_counter] = 0;
		    tmp_counter++;

		}
	    }
	    /*
	     * Put the Ref Neuclotide and Alternative Neuclotide in the dictionary.
	     */
	    int nt_num = temp_content[4].split(",").length;
	    String[] Dictionary = new String[nt_num + 1];
	    Dictionary[0] = temp_content[3];
	    for (int k = 1; k < nt_num + 1; k++) {
		Dictionary[k] = temp_content[4].split(",")[k - 1];
	    }

	    for (int j = 0; j < snp_haplotypes.length; j++) {
		snp_haplotypes[j].AddSNP(temp_content[3], Dictionary[snp_map[j]], this.snp_positions[i]);
	    }

	}

	/*
	 * Statistics for haplotypes
	 */
	HashMap<String, Integer> hap_collection = new HashMap<String, Integer>();

	/*
	 * Here the number of the haplotypes is reduced if this is in test mode.
	 */
	Haplotype[] final_snp_haplotypes;
	if (this.testDatabase) {
	    final_snp_haplotypes = ReduceRandomPortionOfHaplotypes(snp_haplotypes);
	    // final_snp_haplotypes = snp_haplotypes;
	} else {
	    final_snp_haplotypes = snp_haplotypes;
	}

	for (int i = 0; i < final_snp_haplotypes.length; i++) {
	    // if (i == 1214) {
	    // System.out.println("test");
	    // }
	    final_snp_haplotypes[i].GenerateNewHaplotype();
	    if (hap_collection.containsKey(final_snp_haplotypes[i].getHap_seq())) {
		hap_collection.put(final_snp_haplotypes[i].getHap_seq(),
			hap_collection.get(final_snp_haplotypes[i].getHap_seq()) + 1);
	    } else {
		hap_collection.put(final_snp_haplotypes[i].getHap_seq(), 1);
	    }
	}

	Iterator<String> it = hap_collection.keySet().iterator();
	String main_allele_or_not = "";
	int counter = 1;
	while (it.hasNext()) {
	    String key = it.next();

	    if (key.equals(this.ref_seq)) {
		main_allele_or_not = " ***Main allele***";
	    }

	    System.out.println(
		    "Allele" + counter + " " + key + " count is " + hap_collection.get(key) + main_allele_or_not);
	    this.final_result += "Allele" + counter + "p" + " " + key + " count is " + hap_collection.get(key)
		    + main_allele_or_not + System.lineSeparator();
	    main_allele_or_not = "";
	    counter++;
	}
	if (this.sample_map.size() != 0) {
	    this.printHaplotypes(final_snp_haplotypes);
	}
	System.out.println("end");
    }

    public String getFinal_result() {
	return this.final_result;
    }

    private Haplotype[] ReduceRandomPortionOfHaplotypes(Haplotype[] raw_data) {
	Random rand = new Random();
	int randomNum = rand.nextInt(5007);//generate integers range 0 to 5006
	ArrayList<Integer> random_position_to_remove = new ArrayList<Integer>();
	int size = (int) (5008 * this.percentage_of_alleles_to_reduce);
	while (random_position_to_remove.size() != size) {
	    rand = new Random();
	    randomNum = rand.nextInt(5007);
	    if (!random_position_to_remove.contains(randomNum)) {
		random_position_to_remove.add(randomNum);
	    }
	}

	for (int i = 0; i < random_position_to_remove.size(); i++) {
	    int position = random_position_to_remove.get(i);
	    raw_data[position] = null;
	}

	Haplotype[] results = new Haplotype[5008 - random_position_to_remove.size()];
	int j = 0;
	for (int i = 0; i < raw_data.length; i++) {
	    if (raw_data[i] != null) {
		results[j] = raw_data[i];
		j++;
	    }

	}

	return results;
    }

    /*
     * Write the haplotyes of each individual to folder.
     */
    private void printHaplotypes(Haplotype[] haps) throws IOException {
	FileWriter fw = new FileWriter(new File(this.folder_path + "/" + this.genename + "_alleles_info"));
	boolean reversed = false;
	if (Integer.parseInt(rev_info) == 1) {
	    reversed = false;
	} else {
	    reversed = true;
	}
	for (int i = 0; i < haps.length; i++) {
	    Haplotype hap = haps[i];
	    String one_line = "";
	    String seq = hap.getAllele();
	    String rev_seq = "";
	    if (reversed) {
		if (seq == null) {
		    System.out.println("there may be an allele missing!!!!!!!");
		    seq = hap.getHap_seq();
		}
		rev_seq = new ReverseCompliment(seq).getReverseCompliment();
	    } else {
		rev_seq = seq;
	    }
	    if (i % 2 == 0) {
		one_line += hap.getID() + "_chrA" + "\t" + rev_seq + "\t" + hap.getOrigin() + System.lineSeparator();

	    } else {
		one_line += hap.getID() + "_chrB" + "\t" + rev_seq + "\t" + hap.getOrigin() + System.lineSeparator();
	    }
	    fw.write(one_line);
	}
	fw.close();
    }
}
/*
Reverse Compliment:
String input: GGGGaaaaaaaatttatatat
String output: atatataaattttttttCCCC

1. Now, convert each character from input to another character according to the following mapping
public Character execute(Character a) {
'a' -> 't'
'A '-> 'T'
't' -> 'a'
'T' -> 'A'
'g' -> 'c'
'G' -> 'C'
'c' -> 'g'
'C' -> 'G'
'any other' -> 'N' which means i think undefined

if we apply to input then, "CCCCttttttttaaatatata"

2. Revese the result for final output i.e. atatataaattttttttCCCC

references
https://github.com/gifford-lab/GEM/blob/master/src/edu/mit/csail/cgs/ewok/verbs/motifs/ReverseComplement.java
https://www.bioinformatics.org/sms/rev_comp.html
*/
