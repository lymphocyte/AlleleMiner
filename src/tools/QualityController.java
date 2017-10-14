/*
 * The tool is mainly used to trim the alleles with low frequency in the whole number.
 * it will:
 * 1. generate a txt file with the newly discorvered alleles.
 * 2. output a nested hashmap which contains the allels for further use.
 * 
 */
package tools;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.util.ArrayList;
import java.util.HashMap;

public class QualityController {

    private int min_number_of = 1;
    // private ReverseCompliment reverse = new ReverseCompliment();
    // private HashMap<String, HashMap<String, String>> allele_container = new
    // HashMap<String, HashMap<String, String>>();
    ArrayList<String> allele_container = new ArrayList<String>();
    private boolean reversed = false;

    public QualityController(String input_path, String output_path, boolean reversed, int threshold)
	    throws IOException {
	this.reversed = reversed;
	min_number_of = threshold;
	Process(input_path, output_path);

    }

    private void Process(String input_path, String output) throws IOException {
	BufferedReader br = new BufferedReader(new FileReader(input_path));
	FileWriter fw = new FileWriter(new File(output));
	String line = br.readLine();
	String gene_name = "";
	String location = "";
	int allele_counter = 0;
	while (line != null) {
	    /*
	     * If this is the header of the allele.
	     */
	    if (!line.startsWith("Allele")) {
		String[] tmp_header = line.split(" ");
		gene_name = tmp_header[0];
		location = tmp_header[3] + ":" + tmp_header[5];
		allele_counter = 0;
		if (Integer.parseInt(tmp_header[6]) == 1) {
		    this.reversed = false;
		} else {
		    this.reversed = true;
		}
	    } else {
		String[] tmp_allele = line.split(" ");
		/*
		 * Check if the count of the allele exceed the threshold.
		 */
		int count = Integer.parseInt(tmp_allele[4]);
		if (count >= this.min_number_of) {
		    /*
		     * Check if the allele seq contains copy number variation or
		     * not
		     */

		    if (!tmp_allele[1].contains("<")) {

			allele_counter++;
			String allele_fullname = gene_name + "*" + allele_counter;
			String seq = "";
			if (reversed) {
			    ReverseCompliment reverse = new ReverseCompliment(tmp_allele[1]);
			    seq = reverse.getReverseCompliment();
			    // seq = tmp_allele[1];
			} else {
			    seq = tmp_allele[1];
			}
			fw.write(">" + allele_fullname + System.lineSeparator());
			fw.write(seq + System.lineSeparator());
			allele_container.add(">" + allele_fullname + "," + seq + "," + count);
			if (tmp_allele.length == 6) {
			    // allele_fullname_and_count += ",main";
			}

		    } else {
			// continue;
		    }
		}
	    }
	    line = br.readLine();
	}

	fw.close();
	br.close();
    }

    public void SetMinNumber(int num) {
	this.min_number_of = num;
    }

    public ArrayList<String> getAlleles() {
	return this.allele_container;
    }
}
