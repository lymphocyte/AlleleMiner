package input;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileInputStream;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.io.InputStream;
import java.io.InputStreamReader;
import java.io.Reader;
import java.io.UnsupportedEncodingException;
import java.net.MalformedURLException;
import java.net.URL;
import java.util.ArrayList;
import java.util.Properties;
import java.util.zip.GZIPInputStream;

public class ReadVCFfile {
    private int start_index;
    private int end_index;
    private String chr;
    // Properties props;
    String VCF_file_prefix;
    private String result_seq = "";
    private Properties props;
    private String local_vcf_file;
    private StringBuilder result_str_builder = new StringBuilder();
    private ArrayList<String> Indexed_Str = new ArrayList<String>();
    private String[] sample_ids;
    private String hap_path;

    public ReadVCFfile(int start_index, int end_index, String chr, Properties props) throws IOException {
	this.start_index = start_index;
	this.end_index = end_index;
	this.chr = chr;
	this.props = props;
	VCF_file_prefix = props.getProperty("vdj_path");
	this.local_vcf_file = props.getProperty("vcf_path");

	if (props.getProperty("vcf_source").contains("download")) {
	    DownloadFile();
	} else if (props.getProperty("vcf_source").contains("local")) {
	    LoadFromLocalFile();
	} else {
	    LoadFromTextFile();
	}

    }

    private void DownloadFile() throws IOException {
	String url = VCF_file_prefix.split(",")[0] + "chr" + this.chr + VCF_file_prefix.split(",")[1];
	URL link = null;// = new URL("");

	try {
	    link = new URL(url);
	} catch (MalformedURLException e) {
	    // TODO Auto-generated catch block
	    System.out.println("The vdj file doesn't exist in the 1000 genome database");
	    System.exit(0);
	}

	InputStream gzipStream = new GZIPInputStream(link.openStream());

	try {
	    gzipStream = new GZIPInputStream(link.openStream());
	} catch (IOException e) {
	    // TODO Auto-generated catch block
	    System.out.println("The vdj file doesn't exist in the 1000 genome database");
	    System.exit(0);
	}

	Reader decoder = new InputStreamReader(gzipStream, "UTF-8");
	BufferedReader in = new BufferedReader(decoder);

	/*
	 * Read the content in that URL to a string
	 */
	String text = in.readLine();
	System.out.println("####Connected to 1000 genome database####");
	System.out.println("####Searching SNPs of the gene####");

	boolean continueornot = true;
	while (text != null && continueornot) {
	    continueornot = FilteringSNP(text);
	    text = in.readLine();
	}
	// out.close();
	in.close();
	this.result_seq = result_str_builder.toString();
	/*
	 * Write down just the sequences with our interesting sequences.
	 */
	// System.out.println(this.result_seq);
	//
	// FileWriter fw = new FileWriter("IGH.txt");
	// fw.write(this.result_seq);
    }

    private void LoadFromLocalFile() throws IOException {
	System.out.println("####Loading local VCF database####");

	BufferedReader in;

	if (local_vcf_file.contains(".gz")) {
	    InputStream fileStream = new FileInputStream(this.local_vcf_file);
	    InputStream gzipStream = new GZIPInputStream(fileStream);

	    Reader decoder = new InputStreamReader(gzipStream, "UTF-8");
	    in = new BufferedReader(decoder);
	} else {
	    in = new BufferedReader(new FileReader(this.local_vcf_file));
	}

	/*
	 * Read the content in that URL to a string
	 */
	String text = in.readLine();
	boolean continueornot = true;
	while (text != null && continueornot) {
	    continueornot = FilteringSNP(text);
	    text = in.readLine();
	}
	// decoder.close();
	in.close();
	// fileStream.close();
	// gzipStream.close();
	// System.out.println(this.result_seq);
	// System.out.println(this.result_seq);
	this.result_seq = result_str_builder.toString();

	// FileWriter fw = new FileWriter("IGH.txt");
	// fw.write(this.result_seq);
	// fw.close();
	//
    }

    private boolean FilteringSNP(String text) {
	boolean signal = true;
	StringBuilder tmp_builder = new StringBuilder();

	/*
	 * We can get the sample ID from this line.
	 */
	if (text.contains("#CHROM")) {
	    this.sample_ids = text.split("FORMAT")[1].substring(1, text.split("FORMAT")[1].length() - 1).split("\t");
	    this.sample_ids[this.sample_ids.length - 1] = this.sample_ids[this.sample_ids.length - 1]
		    + text.charAt(text.length() - 1);
	}
	if (!text.contains("#")) {
	    int location = Integer.parseInt(text.split("\t")[1]);
	    // System.out.println(location);

	    /*
	     * If there is only one alternative nt one the position.
	     */
	    double AlleFrequency = 0;
	    if (!text.split("\t")[4].contains(",")) {
		// check what we have here.
		String tmp_test = text.split("\t")[7].split(";")[1].split("=")[1];
		if (!tmp_test.contains("0|0")) {
		    AlleFrequency = Double.parseDouble(text.split("\t")[7].split(";")[1].split("=")[1]);
		} else {
		    // System.out.println(text);

		}

		// System.out.println(location);

	    } else {
		String[] snps = text.split("\t")[4].split(",");
		String[] frequencies = text.split("\t")[7].split(";")[1].split("=")[1].split(",");
		int number_of_snps = snps.length;
		for (int i = 0; i < number_of_snps; i++) {
		    if (Double.parseDouble(frequencies[i]) > AlleFrequency) {
			AlleFrequency = Double.parseDouble(frequencies[i]);
		    }
		}
	    }

	    if (location >= this.start_index
		    && AlleFrequency > Double.parseDouble(this.props.getProperty("allele_threshold"))) {
		// System.out.println(text);
		if (location <= this.end_index) {

		    // System.out.println(text);
		    // result_seq += text + System.lineSeparator();
		    tmp_builder.append(text + System.lineSeparator());
		    Indexed_Str.add(text);
		}
	    }

	    if (location > this.end_index) {
		signal = false;
	    }
	}

	return signal;
    }

    public String getRawResult() {
	return this.result_seq;
    }

    public ArrayList<String> getIndexedResults() {
	return this.Indexed_Str;
    }

    /*
     * Just for test, using txt file instead of vcf file. text file contains the
     * snps within interesting regions.
     * 
     */
    private void LoadFromTextFile() throws IOException {
	System.out.println("####Loading local VCF database####");
	InputStream fileStream = new FileInputStream(props.getProperty("filtered_vcf_path"));

	Reader decoder = new InputStreamReader(fileStream);
	BufferedReader in = new BufferedReader(decoder);

	/*
	 * Read the content in that URL to a string
	 */
	String text = in.readLine();
	boolean continueornot = true;
	// System.out.println("SNPs between
	// "+this.start_index+":"+this.end_index);

	while (text != null && continueornot) {
	    continueornot = FilteringSNP(text);
	    text = in.readLine();
	}
	// out.close();
	fileStream.close();
	decoder.close();
	in.close();
	this.result_seq = result_str_builder.toString();

	// System.out.println(this.result_seq);
	// System.out.println(this.result_seq);

    }

    public String[] getSampleIDs() {
	if (this.sample_ids != null) {
	    return this.sample_ids;
	} else {
	    System.out.println("Header file missing");
	    System.exit(0);
	    return null;
	}
    }
}
