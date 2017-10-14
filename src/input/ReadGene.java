package input;

import java.io.IOException;
import java.sql.ResultSet;
import java.sql.SQLException;
import java.util.Properties;

import javax.xml.parsers.ParserConfigurationException;

import org.xml.sax.SAXException;

import tools.DAO;

public class ReadGene {
    private String input;
    private String input_type;
    private String gene_name;
    private String chr;
    private int start_index;
    private int end_index;
    private String DNA_seq;
    Properties prop;
    private String reverse_info = null;

    /*
     * Input type can be genename, or location. If Input is location, it should
     * be like gene_name,chr,start_pos,end_pos
     */
    public ReadGene(String input, String inpuType, Properties prop, String reverse_info) throws SQLException {
	this.input = input;
	this.input_type = inpuType;
	this.prop = prop;

	this.reverse_info = reverse_info;

	process();
    }

    public String get_reverse_info() {
	return this.reverse_info;
    }

    private void process() throws SQLException {
	String[] tmp_location = new String[3];
	if (this.input_type == "genename") {
	    gene_name = this.input;
	    tmp_location = retrieveLocationInfoFromGeneName().split(",");
	    if (tmp_location[0].contains("_")) {
		tmp_location[0] = tmp_location[0].split("_")[0];
	    } else {
	    }
	} else {
	    tmp_location = input.split(",");
	}

	/*
	 * After we get the location information of such gene, we extract the
	 * gene sequence.
	 */
	String seq = "";
	try {
	    GeneLocationExtractor gle = new GeneLocationExtractor(tmp_location[0], Integer.parseInt(tmp_location[1]),
		    Integer.parseInt(tmp_location[2]), this.prop);
	    seq = gle.GetDNAseq();
	} catch (NumberFormatException | ParserConfigurationException | SAXException | IOException e) {
	    // TODO Auto-generated catch block
	    e.printStackTrace();
	}

	this.chr = tmp_location[0];
	this.start_index = Integer.parseInt(tmp_location[1]);
	this.end_index = Integer.parseInt(tmp_location[2]);
	this.DNA_seq = seq;
    }

    private String retrieveLocationInfoFromGeneName() throws SQLException {
	DAO dao = new DAO(prop.getProperty("mysql_host"), prop.getProperty("mysql_user_name"), "");
	String query = "select name,chrom,txStart,txEnd from refGene where name2 = '" + gene_name + "'";
	ResultSet result = dao.search(query);
	System.out.println("####Connected to the reference gene database####");
	if (result == null) {
	    System.out.println(
		    "There is no such gene name, please retype gene name or input the location of this gene.(Be aware the reference genome we use here is chr37)");
	    System.exit(0);
	}

	String result_seq = "";
	while (result.next()) {
	    result_seq = result.getString(2) + "," + result.getInt(3) + "," + result.getInt(4);
	}
	// System.out.println(result_seq);
	if (result_seq == "") {
	    System.out.println(
		    "Couldn't find this gene, try retype your gene name, or type in the location of the gene.");
	    System.exit(0);
	}

	return result_seq;
    }

    public int GetStartIndex() {
	return this.start_index;
    }

    public int GetEndIndex() {
	return this.end_index;
    }

    public String GetDNASeq() {
	return this.DNA_seq;
    }

    public String GetChr() {
	return this.chr;
    }
}
