/// Class to compare seqs with IMGT.
package algorithm; //file is in package algorithm

//including built-in classes
import java.io.BufferedReader;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashMap;
import java.util.Iterator;

//including some other package classes
import org.biojava.nbio.core.exceptions.CompoundNotFoundException;
import org.biojava.nbio.core.sequence.DNASequence;
import org.biojava.nbio.core.sequence.template.Sequence;
import org.biojava.nbio.alignment.Alignments;
import org.biojava.nbio.alignment.Alignments.PairwiseSequenceAlignerType;
import org.biojava.nbio.alignment.SimpleGapPenalty;
import org.biojava.nbio.alignment.SimpleSubstitutionMatrix;
import org.biojava.nbio.alignment.SubstitutionMatrixHelper;
import org.biojava.nbio.alignment.template.SequencePair;
import org.biojava.nbio.alignment.template.SubstitutionMatrix;
import org.biojava.nbio.core.sequence.ProteinSequence;
import org.biojava.nbio.core.sequence.compound.AminoAcidCompound;
import org.biojava.nbio.core.sequence.compound.NucleotideCompound;
import org.biojava.nbio.core.sequence.io.FastaReaderHelper;
import tools.GetRefGene;
import tools.ReverseCompliment;

//main class of the file as it is a public class
public class Compare_With_IMGT {

    ArrayList<String> allele_container = new ArrayList<String>(); //an arraylist with access modifier is default which is only accessible in the current package algorithm
    private String[][] Ref_Gene; //2-d array of String class
    private int novel_size;
    private HashMap<String, HashMap<String, String>> allele_container_imgt; //double hash map, key is String type and value is again HashMap of <String=key, String=value>
    private HashMap<String, HashMap<String, String>> allele_container_novel;
    private String output_prefix = "1000_igh_ref";
    private String gene_type = "";
    private String VJ_type = "";
    private int shiftNUM = 0;
    private String all_files = "";
    private String outputFolder = System.getProperty("user.dir") + "/"; //getProperty of the key "user.dir", which means the return the current working directory
    private String output_all_path = "/Users/yuyaxuan/Desktop/lymDB/Alleles/simulations/simu1.fastq";
    private ArrayList<String> function_genes = new ArrayList<String>();

    public Compare_With_IMGT(ArrayList<String> allele_container, String vjtype, String IMGT_genetype, int number)
	    throws IOException, CompoundNotFoundException {
	/*
	->IOException: exception in the streaming or in simple words file handling types
	->CompoundNotFound: If compound statement does not exist, etc
	*/
	this.allele_container = allele_container;
	output_all_path = "/Users/yuyaxuan/Desktop/lymDB/Alleles/simulations/simu" + number + ".fastq";
	this.VJ_type = vjtype;
	GetRefGene getRef = new GetRefGene("/Users/yuyaxuan/Documents/workspace/AlleleMiner/src/input");
	ArrayList tmpAry = new ArrayList();
	function_genes = getGenenames("/Users/yuyaxuan/Desktop/lymDB/Alleles/germline_TRAV.txt");
	this.gene_type = IMGT_genetype;

	switch (IMGT_genetype) {
	    case "IGH":
		// System.out.println("there are " +
		// getRef.get_Hm_IGHV().length);
		if (VJ_type.contains("V")) {
		    tmpAry.add(getRef.get_Hm_IGHV()); //add data in the arrayList tmpAry
		}
		if (VJ_type.contains("J")) {
		    tmpAry.add(getRef.get_Hm_IGHJ());
		}
		if (VJ_type.contains("D")) {
		    tmpAry.add(getRef.get_Hm_IGHD());
		}
		// tmpAry.add(getRef.get_Hm_IGHD());
		break;
	    case "IGL":
		if (VJ_type.contains("V")) {
		    tmpAry.add(getRef.get_Hm_IGLV());
		}
		if (VJ_type.contains("J")) {
		    tmpAry.add(getRef.get_Hm_IGLJ());
		}
		break;
	    case "TRB":
		if (VJ_type.contains("V")) {
		    tmpAry.add(getRef.get_Hm_TRBV());
		}
		if (VJ_type.contains("J")) {
		    tmpAry.add(getRef.get_Hm_TRBJ());
		}
		if (VJ_type.contains("D")) {
		    tmpAry.add(getRef.get_Hm_TRBD());
		}
		break;
	    case "TRA":
		if (VJ_type.contains("V")) {
		    tmpAry.add(getRef.get_Hm_TRAV());
		}
		if (VJ_type.contains("J")) {
		    tmpAry.add(getRef.get_Hm_TRAJ());
		}
		break;
	    default:
		break;
	}
	this.Ref_Gene = MergeArrayIntoOne(tmpAry);
	this.MappingGeneData();
	CompareBettwenIMGTandNovel();
	WriteExampleRef();
    }

    private String[][] MergeArrayIntoOne(ArrayList<String[][]> al) { //an array list in which each value is of 2-d string array
	int size = 0;
	for (int i = 0; i < al.size(); i++) { //size of array list means number of indexs in array list
	    size += al.get(i).length; //finding the length of index i which is most probably 2
	}
	String[][] refgene = new String[size][2]; //2 d string matrix having column size is 2 and row size is size
	int counter = 0;
	for (int i = 0; i < al.size(); i++) {
	    for (int j = 0; j < al.get(i).length; j++) {
		refgene[counter] = al.get(i)[j]; //find the ith string matrix and the jth row of that matrix
		counter++;
	    }
	}

	return refgene;
    }

    private void CompareBettwenIMGTandNovel() {

    }

	
/*
	 * Following function Convert the novel allele data into a nested hashmap
	 * 
	 * gene - allele - info
	 */	
	
    private void MappingGeneData() throws CompoundNotFoundException, IOException {
	
	HashMap<String, HashMap<String, String>> allele_container = new HashMap<String, HashMap<String, String>>();
	for (int i = 0; i < this.allele_container.size(); i++) {
	    String seq = this.allele_container.get(i);
	    String[] tmp = seq.split(",");
	    String allele_name = tmp[0] + "p";
	    String allele_seq = tmp[1];
	    String genename = allele_name.split("\\*")[0].substring(1); //split allele_name on * and retrieve the 0 index string and the get the subtring from index 1 to rest of the string
	    if (genename.contains(VJ_type)) {
		if (allele_container.containsKey(genename)) { //if genename is present in allele_container then return true otherwise false
		    HashMap<String, String> tmp_inside_allele_container = allele_container.get(genename);
		    tmp_inside_allele_container.put(allele_name, allele_seq); //place allele_name as key and allele_seq as value
		    allele_container.put(genename, tmp_inside_allele_container); //place genename as key and tmp_inside_allele_container as value
		} else {
		    HashMap<String, String> tmp_inside_allele_container = new HashMap<String, String>();
		    tmp_inside_allele_container.put(allele_name, allele_seq);
		    allele_container.put(genename, tmp_inside_allele_container);
		}
	    }
	}
	int counter = 0;
	Iterator<String> it = allele_container.keySet().iterator();
	/*
	->keySet(): return the keys of allele_container in a form of Set
	->iterator(): generate a Interator against the values present in the Set
	*/
	while (it.hasNext()) { //if there is any next element to iterate then return true and move to the first location then second and so on otherwise false
	    String genename = it.next(); //retrieve the next element
	    // System.out.println(genename + ":");
	    Iterator<String> it_inside = allele_container.get(genename).keySet().iterator();
	    while (it_inside.hasNext()) {
		String allelename = it_inside.next();
		// System.out.println("\t" + allelename + ":" +
		// allele_container.get(genename).get(allelename));
		counter++;
	    }
	}

	// this.novel_size = counter;
	/*
	 * Covert IMGT reference gene to HashMap based.
	 */
	HashMap<String, HashMap<String, String>> allele_container_2 = new HashMap<String, HashMap<String, String>>();

	for (int i = 0; i < this.Ref_Gene.length; i++) {

	    String allele_name = this.Ref_Gene[i][0];
	    String allele_seq = this.Ref_Gene[i][1];
	    String genename = this.Ref_Gene[i][0].split("\\*")[0].substring(1);
	    if (genename.contains(VJ_type)) {
		if (allele_container_2.containsKey(genename)) {
		    HashMap<String, String> tmp_inside_allele_container = allele_container_2.get(genename);
		    tmp_inside_allele_container.put(allele_name, allele_seq);
		    allele_container_2.put(genename, tmp_inside_allele_container);
		} else {
		    HashMap<String, String> tmp_inside_allele_container = new HashMap<String, String>();
		    tmp_inside_allele_container.put(allele_name, allele_seq);
		    allele_container_2.put(genename, tmp_inside_allele_container);
		}
	    }
	}

	allele_container_imgt = allele_container_2;
	allele_container_novel = allele_container;
	// System.out.println("There are " + allele_container_2.size() + " genes
	// in IMGT");
	// System.out.println("There are " + allele_container.size() + " genes
	// in novel");

	ComparingIMGTwithNovel();
	// System.out.println(total);
    }

    private void ComparingIMGTwithNovel() throws CompoundNotFoundException, IOException {
	ArrayList<String> novel_results = new ArrayList<String>(); //array list of string type
	int total = 0;
	int sub = 0;
	int another = 0;
	int novel_matched = 0;
	int novel_total = 0;
	int num_no_gene_allele = 0;
	HashMap<String, HashMap<String, String>> allele_both = new HashMap<String, HashMap<String, String>>();
	HashMap<String, HashMap<String, String>> allele_left_from_novel = new HashMap<String, HashMap<String, String>>();
	HashMap<String, HashMap<String, String>> allele_left_from_imgt = new HashMap<String, HashMap<String, String>>();
	/*
	 * IMGT allele information
	 */
	HashMap<String, ArrayList<String>> IMGT_allele_result_info = new HashMap<String, ArrayList<String>>();
	String result_imgt_not_in_novel = "";
	String result_imgt_and_in_novel = "";
	String resut_imgt_ale_not_found = "";
	String resut_imgt_gene_not_found = "";
	ArrayList<String> gene_names = new ArrayList<String>();
	ArrayList<String> ale_names = new ArrayList<String>();

	/*
	 * Starting from each IMGT sequence to check which IMGT sequences are
	 * not in novel sequences.
	 */
	Iterator<String> it = allele_container_imgt.keySet().iterator();
	int counter_gene_num = 0;
	int non_existed_imgt_gene_num = 0;
	int counterofrealalleles = 0;
	while (it.hasNext()) {
	    counter_gene_num++;
	    String genename = it.next();
	    // System.out.println(genename + ":");
	    /*
	     * Novel genes inferd from 1000 genome data may include different
	     * names, so we first check if the gene name from IMGT does exist in
	     * novel gene sets.
	     */
	    if (this.allele_container_novel.containsKey(genename)) {

		HashMap<String, String> tmp = this.allele_container_imgt.get(genename); //return value against key named genename from allele_container_imgt hash map
		Iterator<String> it_inside = tmp.keySet().iterator();
		novel_total += this.allele_container_novel.get(genename).size();
		while (it_inside.hasNext()) {
		    String allelename = it_inside.next();
		    int last_num = 0;
		    if (this.VJ_type.contains("V")) {
			last_num = this.shiftNUM;
		    }
		    String imgt_seq = this.allele_container_imgt.get(genename).get(allelename); /*first return value from allele_container_imgt which is hashmap<String, String> and the
												again get value from it against key names allelename*/  

		    Iterator<String> it_novel = this.allele_container_novel.get(genename).keySet().iterator();
		    int flag = 0;
		    boolean alreadyADD = false;
		    String novel_alleles = "";
		    while (it_novel.hasNext()) {
			String novel_allele_name = it_novel.next();

			String novel_allele = this.allele_container_novel.get(genename).get(novel_allele_name);
			ReverseCompliment reverse = new ReverseCompliment(novel_allele);
			String rev_novel_allele = reverse.getReverseCompliment();
			String imgt_allele = this.allele_container_imgt.get(genename).get(allelename);
			/*
			 * TEST
			 */
			// if (allelename.contains("IGHJ5*02")) {
			// System.out.println("here");
			// }
			/*
			 * END
			 */
			if (novel_allele.contains(imgt_allele.substring(last_num, imgt_allele.length() - last_num))
				|| rev_novel_allele
					.contains(imgt_allele.substring(last_num, imgt_allele.length() - last_num))) {
				/*
				->contains: it returns true if the given string (as a paramter) is present as a substring in the main string (by which contains is called) otherwise return false.
					its method of finding the substring is case-sensitive
				->substring:(int start): it simply retrieves the substring of the main string (by which function is called) from start index to length-1
				->subtring:(int start, int end): it retrieves the substring of the main string (by which function is called) from start index to end-1
				*/
			    if (!alreadyADD) { //boolean variable defined above
				novel_matched++;
				alreadyADD = true;
			    }
			    flag++;
			    if (novel_results.contains(novel_allele_name)) { //find novel_allele_name in novel_results and return true/false accordingly

			    } else {
				novel_results.add(novel_allele_name);
				novel_alleles += novel_allele_name.substring(1) + ",";
			    }
			    // AlignTwoSeq(this.allele_container_novel.get(genename).get(novel_allele_name),
			    // this.allele_container_imgt.get(genename).get(allelename));
			    // both_alleles.put(allelename, value)
			} else {

			}

		    }
		    total++;
		    if (flag == 0) {

			// If we can't find a perfect match of this allele in
			// the novel allele database, then we align them using
			// local alignment algorithm.
			String IMGT_allle_seq = this.allele_container_imgt.get(genename).get(allelename);
			result_imgt_not_in_novel += allelename + System.lineSeparator() + IMGT_allle_seq
				+ System.lineSeparator(); // System.lineSeperator: for new line, in windows \r\n and in UNIX \n
			Iterator<String> aligner_itetrator = this.allele_container_novel.get(genename).keySet()
				.iterator();

			another++;

			/*
			 * check if these alelles are from rearranged or not.
			 */

			if (!this.function_genes.contains(allelename.substring(1))) {
			    // this.function_genes.remove(allelename.substring(1));
			    ArrayList<String> tmp_ary = new ArrayList<String>();
			    tmp_ary.add("This IMGT alelle is retrieved from rearranged sequence");
			    tmp_ary.add("");
			    IMGT_allele_result_info.put(allelename, tmp_ary);
			} else {
			    ArrayList<String> tmp_ary = new ArrayList<String>();
			    tmp_ary.add("This IMGT allele can't be recovered from G1K data");
			    tmp_ary.add("");
			    IMGT_allele_result_info.put(allelename, tmp_ary);
			}
		    } else {
			result_imgt_and_in_novel += allelename + System.lineSeparator()
				+ this.allele_container_imgt.get(genename).get(allelename) + System.lineSeparator();
			sub++;

			/*
			 * Adding this to count the number of real alleles from
			 * IMGT being retrived.
			 */
			ArrayList<String> tmp_ary = new ArrayList<String>();
			tmp_ary.add(" ");
			tmp_ary.add(novel_alleles);
			IMGT_allele_result_info.put(allelename, tmp_ary);
			if (this.function_genes.contains(allelename.substring(1))) {
			    // this.function_genes.remove(allelename.substring(1));
			    // System.out.println(allelename);
			    counterofrealalleles++;
			}
		    }
		}
	    } else {
		HashMap<String, String> tmp = this.allele_container_imgt.get(genename);
		Iterator<String> ale_name = tmp.keySet().iterator();
		while (ale_name.hasNext()) {
		    ArrayList<String> tmp_ary = new ArrayList<String>();
		    tmp_ary.add("This gene is not in the current genomoe build");
		    tmp_ary.add("");
		    IMGT_allele_result_info.put(ale_name.next(), tmp_ary);
		}
		/*
		 * Let's see which genes can't be found in the novel gene set.
		 */
		// System.out.print("\""+genename+"\",");

		/*
		 * We write the alleles that don't belong to any genes to a file
		 * here.
		 */
		non_existed_imgt_gene_num++;
		Iterator<String> it_another = this.allele_container_imgt.get(genename).keySet().iterator();
		while (it_another.hasNext()) {
		    num_no_gene_allele++;
		    String allelename = it_another.next();
		    String genenname = allelename.split("\\*")[0].split(">")[1];
		    resut_imgt_ale_not_found += allelename + System.lineSeparator()
			    + this.allele_container_imgt.get(genename).get(allelename) + System.lineSeparator();
		    if (gene_names.contains(genenname)) {

		    } else {
			gene_names.add(genenname);
			resut_imgt_gene_not_found += genename + "\t";
		    }
		}
		// result_imgt_and_in_novel += allelename +
		// System.lineSeparator()+
		// this.allele_container_imgt.get(genename).get(allelename) +
		// System.lineSeparator();
	    }
	}

	write_result_of_ale_info(IMGT_allele_result_info);
	FileWriter fw = new FileWriter(new File("Alleles_without_gene.txt")); //fileWriter: used for writing file. if file is not present, it will create it. By default it is not in appending mode. Every time a new file will be created
	fw.write(resut_imgt_gene_not_found); //write resut_imgt_gene_not_found in file
	fw.close(); //close the FileWriter object
        //FileWriter and its use can cause IOException, so we have to handle it explicitly
	
	// System.out.println();
	novel_size = novel_total;
	System.out.print(novel_size + ",");
	System.out.println(sub);
	System.out.println(sub + " IMGT alleles found in novel alleles out of " + total + " IMGT alllels");
	System.out.println("The number of alleles that are in IMGT but not stored as the varation of genes in "
		+ num_no_gene_allele);
	System.out.println(
		"There are " + novel_results.size() + " novel alleles out of " + novel_size + " can be found in IMGT!");
	System.out.println("there are " + counterofrealalleles + " out of " + this.function_genes.size());
	// System.out.print(novel_size-novel_results.size()+",");
	// // CompareNovelWithIMGT();
	FingGenesNotinIMGT();

	FileWriter fw_imgt_no = new FileWriter(
		new File(this.outputFolder + this.gene_type + this.VJ_type + "_imgt_not_in_novel.txt"));
	FileWriter fw_imgt_in = new FileWriter(
		new File(this.outputFolder + this.gene_type + this.VJ_type + "_imgt_and_in_novel.txt"));

	fw_imgt_no.write(result_imgt_not_in_novel);
	fw_imgt_in.write(result_imgt_and_in_novel);
	fw_imgt_no.close();
	fw_imgt_in.close();
	System.out.println("There are " + counter_gene_num + " genes in total.");
	System.out.println("There are " + non_existed_imgt_gene_num + " genes can't be found!");
    }

    /*
     * This will filter out the genes(mainly open reading frame and
     * pseudogenes). And stores the
     */
    private void FingGenesNotinIMGT() {
	int counter = 0;
	StringBuilder sb_V = new StringBuilder();
	StringBuilder sb_D = new StringBuilder();
	StringBuilder sb_J = new StringBuilder();

	Iterator<String> it = this.allele_container_novel.keySet().iterator();
	while (it.hasNext()) {
	    String gene = it.next();
	    Iterator<String> it_inside = this.allele_container_imgt.keySet().iterator();
	    if (allele_container_imgt.containsKey(gene)) {
		Iterator<String> it_imgt = allele_container_novel.get(gene).keySet().iterator();
		while (it_imgt.hasNext()) {
		    String allelename = it_imgt.next();
		    String seq = allele_container_novel.get(gene).get(allelename);
		    if (allelename.contains("V")) {
			sb_V.append(allelename); // add allelename at the end of sb_V string
			sb_V.append(System.lineSeparator());
			sb_V.append(seq);
			sb_V.append(System.lineSeparator());
		    } else if (allelename.contains("D")) {
			sb_D.append(allelename);
			sb_D.append(System.lineSeparator());
			sb_D.append(seq);
			sb_D.append(System.lineSeparator());
		    } else if (allelename.contains("J")) {
			sb_J.append(allelename);
			sb_J.append(System.lineSeparator());
			sb_J.append(seq);
			sb_J.append(System.lineSeparator());
		    } else {
			System.out.println("There is something wrong here");
		    }

		}
	    } else {
		counter++;
		// System.out.println(gene);
	    }
	}

    }

    // Compare the novel alleles with IMGT alleles
    private void CompareNovelWithIMGT() throws IOException {

	Iterator<String> it_novel = this.allele_container_novel.keySet().iterator();
	int novel_counts = 0;
	int num_of_non_gene = 0;
	String writeto_file = "";
	int allele_counter = 0;
	int real_novel_size = 0;
	while (it_novel.hasNext()) {
	    String genename_novel = it_novel.next();

	    if (this.allele_container_imgt.containsKey(genename_novel)) {
		Iterator<String> it_novel_allele = this.allele_container_novel.get(genename_novel).keySet().iterator();
		while (it_novel_allele.hasNext()) {
		    allele_counter++;
		    // real_novel_size;
		    String allelename_novel = it_novel_allele.next();
		    int last_num = 0;
		    if (allelename_novel.contains(this.gene_type + "V")) {
			last_num = shiftNUM;
		    }

		    String novel_allele = this.allele_container_novel.get(genename_novel).get(allelename_novel);
		    // ReverseCompliment reverse = new
		    // ReverseCompliment(novel_allele);
		    // String rev_novel_allele_seq =
		    // reverse.getReverseCompliment();
		    // Now we loop IMGT alleles.
		    Iterator<String> it_imgt = this.allele_container_imgt.get(genename_novel).keySet().iterator();
		    boolean flag = false;
		    while (it_imgt.hasNext()) {
			String imgt_allele_name = it_imgt.next();
			String imgt_allele = this.allele_container_imgt.get(genename_novel).get(imgt_allele_name);
			ReverseCompliment reverse = new ReverseCompliment(novel_allele);
			String rev_novel_allele = reverse.getReverseCompliment();
			if (novel_allele.contains(imgt_allele.substring(last_num, imgt_allele.length() - last_num))
				|| rev_novel_allele
					.contains(imgt_allele.substring(last_num, imgt_allele.length() - last_num))
				|| imgt_allele
					.contains(novel_allele.substring(last_num, novel_allele.length() - last_num))
				|| imgt_allele.contains(rev_novel_allele.substring(last_num, rev_novel_allele.length() - last_num))) {
			    // novel_counts++;
			    flag = true;
			}
		    }

		    if (flag) {
			novel_counts++;
		    } else {
			writeto_file += allelename_novel + System.lineSeparator() + novel_allele
				+ System.lineSeparator();

			// writeto_file += allelename_novel + "_rev" +
			// System.lineSeparator() + rev_novel_allele_seq
			// + System.lineSeparator();
		    }
		}
	    } else {
		// System.out.println(genename_novel);
		Iterator it = this.allele_container_novel.get(genename_novel).keySet().iterator();
		while (it.hasNext()) {
		    num_of_non_gene++;
		    it.next();
		}
	    }
	}
	System.out.println("number of alleles that are in novel but not stored as variation of gene in IMGT is "
		+ num_of_non_gene);
	System.out.println(novel_counts + " novel alleles found in IMGT out of " + novel_size);
	FileWriter fw = new FileWriter(
		new File(this.outputFolder + this.gene_type + this.VJ_type + "_novel_not_in_imgt.txt"));
	fw.write(writeto_file);
	all_files = writeto_file;
	fw.close();
    }

    // Using BioJava to do local alignment for two seqs.

    private SequencePair<DNASequence, NucleotideCompound> AlignTwoSeq(String seq1, String seq2)
	    throws CompoundNotFoundException {

	DNASequence dna_seq_1 = new DNASequence(seq1);
	DNASequence dna_seq_2 = new DNASequence(seq2);
	SubstitutionMatrix<NucleotideCompound> matrix = SubstitutionMatrixHelper.getNuc4_4();

	SimpleGapPenalty gapP = new SimpleGapPenalty();
	gapP.setOpenPenalty((short) 5);
	gapP.setExtensionPenalty((short) 2);

	SequencePair<DNASequence, NucleotideCompound> psa = Alignments.getPairwiseAlignment(dna_seq_1, dna_seq_2,
		PairwiseSequenceAlignerType.LOCAL, gapP, matrix);

	return psa;
	// System.out.println(psa);
	// System.out.println(psa.getNumSimilars());
	// System.out.println(psa.getNumIdenticals());
	// System.out.println(psa.getLength());
	// System.out.println(seq1.length());
	// System.out.println(seq2.length());

    }

    /*
     * We write a combination of IMGT and novel database.
     */
    private void WriteExampleRef() throws IOException {
	String all_path = "/Users/yuyaxuan/Desktop/lymDB/Alleles/All_IMGT+New_alleles.txt";
	String write_to_file = "";
	Iterator<String> it = allele_container_imgt.keySet().iterator();
	while (it.hasNext()) {
	    String gene_name = it.next();
	    Iterator<String> it_ale = allele_container_imgt.get(gene_name).keySet().iterator();
	    while (it_ale.hasNext()) {
		String ale_name = it_ale.next();
		write_to_file += ale_name + System.lineSeparator()
			+ this.allele_container_imgt.get(gene_name).get(ale_name) + System.lineSeparator();
	    }
	}

	it = this.allele_container_novel.keySet().iterator();
	while (it.hasNext()) {
	    String gene_name = it.next();
	    Iterator<String> it_ale = allele_container_novel.get(gene_name).keySet().iterator();
	    while (it_ale.hasNext()) {
		String ale_name = it_ale.next();
		write_to_file += ale_name + System.lineSeparator()
			+ this.allele_container_novel.get(gene_name).get(ale_name) + System.lineSeparator();
	    }
	}

	FileWriter fw = new FileWriter(all_path); //all_path is defined above
	fw.write(write_to_file);
	fw.close();
    }

    private ArrayList<String> getGenenames(String filename) throws IOException {
	BufferedReader br = new BufferedReader(new FileReader(filename)); //BufferedReader and FileReader are the object for file reading
	String line = br.readLine();//read a single line upto to \n
	String[] ary = line.split(",");
	return new ArrayList<String>(Arrays.asList(ary));
    }

    private void write_result_of_ale_info(HashMap<String, ArrayList<String>> result) throws IOException {

	String filename_ale_in_current_build = "/Users/yuyaxuan/Desktop/lymDB/ale_result/" + this.gene_type
		+ this.VJ_type + ".csv";
	String filename_ale_not_in_current_build = "/Users/yuyaxuan/Desktop/lymDB/ale_result/" + this.gene_type
		+ "not_in_current_build.csv";

	FileWriter fw = new FileWriter(new File(filename_ale_in_current_build));
	FileWriter fw_not = new FileWriter(new File(filename_ale_not_in_current_build));
	fw.write("IMGT allele" + "\t" + "Allele(s) from Lym1K" + "\t" + "Comments" + System.lineSeparator());
	Iterator<String> it = result.keySet().iterator();
	while (it.hasNext()) {
	    String alename = it.next();
	    if (result.get(alename).get(1).length() > 0) {
		fw.write(alename.substring(1) + "\t"
			+ result.get(alename).get(1).substring(0, result.get(alename).get(1).length() - 1) + "\t"
			+ result.get(alename).get(0) + System.lineSeparator());
	    } else {
		// if (result.get(alename).get(0).contains("current")) {
		// fw_not.write(alename.substring(1) + "\t" +
		// result.get(alename).get(1) + "\t"
		// + result.get(alename).get(0) + System.lineSeparator());
		// } else {
		fw.write(alename.substring(1) + "\t" + result.get(alename).get(1) + "\t" + result.get(alename).get(0)
			+ System.lineSeparator());
		// }
	    }
	}

	fw.close();
	fw_not.close();
    }
}
