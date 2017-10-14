
package tools;

import java.io.BufferedReader;
import java.io.FileInputStream;
import java.io.IOException;
import java.io.InputStream;
import java.io.InputStreamReader;
import java.net.URL;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.Iterator;

public class GetRefGene {

    private String[][] Hm_IGHV;
    private String[][] Hm_TRBV_NoPoly;
    private String[][] Hm_IGHV_NO_POLY;
    private String[][] Hm_IGHJ_NO_POLY;
    private String[][] Hm_IGHJ;
    private String[][] Hm_IGHD;
    private String[][] Hm_IGKV;
    private String[][] Hm_IGKJ;
    private String[][] Hm_IGLV;
    private String[][] Hm_IGLJ;
    private String[][] Ms_IGHV;
    private String[][] Ms_IGHJ;
    private String[][] Ms_IGHD;
    private String[][] Ms_IGKV;
    private String[][] Ms_IGKJ;
    private String[][] Ms_IGLV;
    private String[][] Ms_IGLJ;
    private String[][] Ms_IGJ;
    private String[][] Ms_TRBV;
    private String[][] Ms_TRBJ;
    private String[][] Ms_TRBD;
    private String[][] Ms_TRAV;
    private String[][] Ms_TRAJ;
    private String[][] Ms_TRDD;
    private String[][] Ms_TRDV;
    private String[][] Ms_TRDJ;
    private String[][] Ms_TRGV;
    private String[][] Ms_TRGJ;
    private String[][] Hm_IGV;
    private String[][] Hm_IGJ;
    private String[][] Hm_TRBV;
    private String[][] Hm_TRBD;
    private String[][] Hm_TRBJ;
    private String[][] Hm_TRAV;
    private String[][] Hm_TRAJ;
    private String[][] Hm_TRGV;
    private String[][] Hm_TRGJ;
    private String[][] Hm_TRDV;
    private String[][] Hm_TRDJ;
    private String[][] Hm_TRDD;
    private String[][] Hm_TRBJ_NO_poly;
    private String Ref_Input_Folder = "";

    /*
     * Constructor without
     */

    public GetRefGene(String seq) {
	if (seq.equals("")) {
	    this.Ref_Input_Folder = seq;
	} else {
	    if (!(seq.charAt(seq.length() - 1) + "").equals("/")) {
		seq = seq + "/";
	    }
	    this.Ref_Input_Folder = seq;
	}
    }

    public String[][] get_Ms_TRBV() throws IOException {
	Ms_TRBV = ConvertHashMapToArrary(load(this.Ref_Input_Folder + "mouse_TRBV_region.fasta"));
	return this.Ms_TRBV;
    }

    public String[][] get_Ms_TRBJ() throws IOException {
	Ms_TRBJ = ConvertHashMapToArrary(load(this.Ref_Input_Folder + "mouse_TRBJ_region.fasta"));
	return this.Ms_TRBJ;
    }

    public String[][] get_Ms_TRBD() throws IOException {
	Ms_TRBD = ConvertHashMapToArrary(load(this.Ref_Input_Folder + "mouse_TRBD"));
	return this.Ms_TRBD;
    }

    public String[][] get_Ms_TRAV() throws IOException {
	Ms_TRAV = ConvertHashMapToArrary(load(this.Ref_Input_Folder + "mouse_TRAV_region.fasta"));
	return this.Ms_TRAV;
    }

    public String[][] get_Ms_TRAJ() throws IOException {
	Ms_TRAJ = ConvertHashMapToArrary(load(this.Ref_Input_Folder + "mouse_TRBJ_region.fasta"));
	return this.Ms_TRAJ;
    }

    public String[][] get_Ms_TRDD() throws IOException {
	Ms_TRDD = ConvertHashMapToArrary(load(this.Ref_Input_Folder + "mouse_TRDD"));
	return this.Ms_TRDD;
    }

    public String[][] get_Ms_TRDV() throws IOException {
	Ms_TRDV = ConvertHashMapToArrary(load(this.Ref_Input_Folder + "mouse_TRDV_region.fasta"));
	return this.Ms_TRDV;
    }

    public String[][] get_Ms_TRDJ() throws IOException {
	Ms_TRDJ = ConvertHashMapToArrary(load(this.Ref_Input_Folder + "mouse_TRDJ_region.fasta"));
	return this.Ms_TRDJ;
    }

    public String[][] get_Ms_TRGV() throws IOException {
	Ms_TRGV = ConvertHashMapToArrary(load(this.Ref_Input_Folder + "mouse_TRGV_region.fasta"));
	return this.Ms_TRGV;
    }

    public String[][] get_Ms_TRGJ() throws IOException {
	Ms_TRGJ = ConvertHashMapToArrary(load(this.Ref_Input_Folder + "mouse_TRGJ_region.fasta"));
	return this.Ms_TRGJ;
    }

    public String[][] get_Hm_TRGV() throws IOException {
	Hm_TRGV = ConvertHashMapToArrary(load(this.Ref_Input_Folder + "human_TRGV_region.fasta"));
	return this.Hm_TRGV;
    }

    public String[][] get_Hm_TRGJ() throws IOException {
	Hm_TRGJ = ConvertHashMapToArrary(load(this.Ref_Input_Folder + "human_TRGJ_region.fasta"));
	return this.Hm_TRGJ;
    }

    public String[][] get_Hm_TRDV() throws IOException {
	Hm_TRDV = ConvertHashMapToArrary(load(this.Ref_Input_Folder + "human_TRDV_region.fasta"));
	return this.Hm_TRDV;
    }

    public String[][] get_Hm_TRDJ() throws IOException {
	Hm_TRDJ = ConvertHashMapToArrary(load(this.Ref_Input_Folder + "human_TRDJ_region.fasta"));
	return this.Hm_TRDJ;
    }

    public String[][] get_Hm_TRBD() throws IOException {
	Hm_TRBD = ConvertHashMapToArrary(load(this.Ref_Input_Folder + "human_TRBD_region.fasta"));
	return this.Hm_TRBD;
    }

    public ArrayList<String> get_Hm_TRBD_genemaes() throws IOException {
	ArrayList<String> results = new ArrayList<String>();
	Hm_TRBD = ConvertHashMapToArrary(load(this.Ref_Input_Folder + "human_TRBD"));
	for (int i = 0; i < Hm_TRBD.length; i++) {
	    if (results.contains(Hm_TRBD[i])) {

	    } else {
		results.add(Hm_TRBD[i][0]);
	    }
	}
	return results;
    }

    public String[][] get_Hm_IGJ() {
	return this.Hm_IGJ;
    }

    public String[][] get_Hm_IGV() {
	return this.Hm_IGV;
    }

    public String[][] get_Ms_IGHV() throws IOException {
	Ms_IGHV = ConvertHashMapToArrary(load(this.Ref_Input_Folder + "mouse_IGHV"));
	return this.Ms_IGHV;
    }

    public String[][] get_Ms_IGHD() throws IOException {
	Ms_IGHD = ConvertHashMapToArrary(load(this.Ref_Input_Folder + "mouse_IGHD"));
	return this.Ms_IGHD;
    }

    public String[][] get_Ms_IGHJ() throws IOException {
	Ms_IGHJ = ConvertHashMapToArrary(load(this.Ref_Input_Folder + "mouse_IGHJ"));
	return this.Ms_IGHJ;
    }

    public String[][] get_Ms_IGKV() {
	return this.Ms_IGKV;
    }

    public String[][] get_Ms_IGKJ() {
	return this.Ms_IGKJ;
    }

    public String[][] get_Hm_TRDD() throws IOException {
	Hm_TRDD = ConvertHashMapToArrary(load(this.Ref_Input_Folder + "human_TRDD"));
	return this.Hm_TRDD;
    }

    public String[][] get_Ms_IGLV() throws IOException {
	Ms_IGLV = ConvertHashMapToArrary(load(this.Ref_Input_Folder + "mouse_IGLV"));
	return this.Ms_IGLV;
    }

    public String[][] get_Ms_IGLJ() throws IOException {
	Ms_IGLJ = ConvertHashMapToArrary(load(this.Ref_Input_Folder + "mouse_IGLJ"));
	return this.Ms_IGLJ;
    }

    public String[][] get_Hm_TRAV() throws IOException {
	Hm_TRAV = ConvertHashMapToArrary(load(this.Ref_Input_Folder + "human_full_TRAV"));
	return this.Hm_TRAV;
    }

    public String[][] get_Hm_TRAJ() throws IOException {
	Hm_TRAJ = ConvertHashMapToArrary(load(this.Ref_Input_Folder + "human_full_TRAJ"));
	return this.Hm_TRAJ;
    }

    public String[][] get_Hm_IGHV() throws IOException {
	Hm_IGHV = ConvertHashMapToArrary(load(this.Ref_Input_Folder + "human_IGHV"));
	return this.Hm_IGHV;
    }

    public ArrayList<String> get_Hm_IGHV_genemaes() throws IOException {
	ArrayList<String> results = new ArrayList<String>();
	Hm_IGHV = ConvertHashMapToArrary(load(this.Ref_Input_Folder + "human_IGHV"));
	for (int i = 0; i < Hm_IGHV.length; i++) {
	    if (results.contains(Hm_IGHV[i])) {

	    } else {
		results.add(Hm_IGHV[i][0]);
	    }
	}
	return results;
    }

    public String[][] get_Hm_IGHV_NO_POLY() throws IOException {
	Hm_IGHV_NO_POLY = ConvertHashMapToArrary(load(this.Ref_Input_Folder + "human_IGHV_NO_POLY.fasta"));
	return this.Hm_IGHV_NO_POLY;
    }

    public String[][] get_Hm_IGHD() throws IOException {
	Hm_IGHD = ConvertHashMapToArrary(load(this.Ref_Input_Folder + "human_IGHD"));
	return this.Hm_IGHD;
    }

    public ArrayList<String> get_Hm_IGHD_genemaes() throws IOException {
	ArrayList<String> results = new ArrayList<String>();
	Hm_IGHD = ConvertHashMapToArrary(load(this.Ref_Input_Folder + "human_IGHD"));
	for (int i = 0; i < Hm_IGHD.length; i++) {
	    if (results.contains(Hm_IGHD[i])) {

	    } else {
		results.add(Hm_IGHD[i][0]);
	    }
	}
	return results;
    }

    public String[][] get_Hm_IGHJ() throws IOException {
	Hm_IGHJ = ConvertHashMapToArrary(load(this.Ref_Input_Folder + "human_IGHJ"));
	return this.Hm_IGHJ;
    }

    public ArrayList<String> get_Hm_IGHJ_genemaes() throws IOException {
	ArrayList<String> results = new ArrayList<String>();
	Hm_IGHJ = ConvertHashMapToArrary(load(this.Ref_Input_Folder + "human_IGHJ"));
	for (int i = 0; i < Hm_IGHJ.length; i++) {
	    if (results.contains(Hm_IGHJ[i])) {

	    } else {
		results.add(Hm_IGHJ[i][0]);
	    }
	}
	return results;
    }

    public String[][] get_Hm_IGHJ_NO_POLY() throws IOException {
	Hm_IGHJ_NO_POLY = ConvertHashMapToArrary(load(this.Ref_Input_Folder + "human_IGJ_NO_POLY"));
	return this.Hm_IGHJ_NO_POLY;
    }

    public String[][] get_Hm_IGKV() {
	return this.Hm_IGKV;
    }

    public String[][] get_Hm_IGKJ() {
	return this.Hm_IGKJ;
    }

    public String[][] get_Hm_IGLV() throws IOException {
	Hm_IGLV = ConvertHashMapToArrary(load(this.Ref_Input_Folder + "human_IGLV"));
	return this.Hm_IGLV;
    }

    public String[][] get_Hm_IGLJ() throws IOException {
	Hm_IGLJ = ConvertHashMapToArrary(load(this.Ref_Input_Folder + "human_IGLJ"));
	return this.Hm_IGLJ;
    }

    public String[][] get_Hm_full_TRBV() throws IOException {
	Hm_TRBV = ConvertHashMapToArrary(load(this.Ref_Input_Folder + "human_full_TRBV"));
	return this.Hm_TRBV;
    }

    public String[][] get_Hm_TRBV() throws IOException {
	Hm_TRBV = ConvertHashMapToArrary(load(this.Ref_Input_Folder + "human_TRBV_region.fasta"));
	// String seq =
	// load(this.Ref_Input_Folder+"human_TRBV_region.fasta").get(">TRBV11-2*01");
	return this.Hm_TRBV;
    }

    public ArrayList<String> get_Hm_TRBV_genemaes() throws IOException {
	ArrayList<String> results = new ArrayList<String>();
	Hm_TRBV = ConvertHashMapToArrary(load(this.Ref_Input_Folder + "human_TRBV"));
	for (int i = 0; i < Hm_TRBV.length; i++) {
	    if (results.contains(Hm_TRBV[i])) {

	    } else {
		results.add(Hm_TRBV[i][0]);
	    }
	}
	return results;
    }

    public String[][] get_Hm_TRBV_NoPoly() throws IOException {
	this.Hm_TRBV_NoPoly = ConvertHashMapToArrary(load(this.Ref_Input_Folder + "TRBV_No_poly.fasta"));
	return this.Hm_TRBV_NoPoly;
    }

    public String[][] get_Hm_TRBJ() throws IOException {
	Hm_TRBJ = ConvertHashMapToArrary(load(this.Ref_Input_Folder + "human_full_TRBJ"));
	return this.Hm_TRBJ;
    }

    public ArrayList<String> get_Hm_TRBJ_genemaes() throws IOException {
	ArrayList<String> results = new ArrayList<String>();
	Hm_TRBJ = ConvertHashMapToArrary(load(this.Ref_Input_Folder + "human_TRBJ"));
	for (int i = 0; i < Hm_TRBJ.length; i++) {
	    if (results.contains(Hm_TRBJ[i])) {

	    } else {
		results.add(Hm_TRBJ[i][0]);
	    }
	}
	return results;
    }

    public String[][] get_Hm_TRBJ_NO_poly() throws IOException {
	this.Hm_TRBJ_NO_poly = ConvertHashMapToArrary(load(this.Ref_Input_Folder + "human_TRBJ_NO_Poly.fasta"));
	return this.Hm_TRBJ_NO_poly;
    }

    public String[][] get_Ms_IGJ() throws IOException {
	Ms_IGJ = ConvertHashMapToArrary(load(this.Ref_Input_Folder + "mouse_IGJ"));
	return this.Ms_IGJ;
    }

    HashMap<String, String> load(String filename) throws IOException {
	InputStream stream;
	if (filename.contains("/")) {
	    stream = new FileInputStream(filename);
	} else {
	    URL urlToDictionary = this.getClass().getResource(filename);
	    stream = urlToDictionary.openStream();
	}

	HashMap<String, String> genes = new HashMap<String, String>();
	BufferedReader br = new BufferedReader(new InputStreamReader(stream));
	// BufferedReader tmp_br =br;
	String strLine;

	String tmp_gene_name = "";

	String seqFile = "";
	while ((strLine = br.readLine()) != null) {
	    seqFile += strLine + System.getProperty("line.separator");
	}

	String[] tmp_seqFile = seqFile.split(System.getProperty("line.separator"));
	boolean tag = true;
	for (int i = 0; i < tmp_seqFile.length; i++) {
	    if (tmp_seqFile[i].startsWith(">")) {
		String name = tmp_seqFile[i];
		int j = i + 1;
		String tmp_seq = "";
		int size = tmp_seqFile.length;
		if (j < size) {
		    while (!tmp_seqFile[j].startsWith(">")) {
			tmp_seq += tmp_seqFile[j];
			j++;
			if (j >= size) {
			    break;
			}
		    }
		}

		genes.put(name, tmp_seq);
		i = j - 1;
	    }

	}

	return genes;

    }

    private String[][] ConvertHashMapToArrary(HashMap<String, String> genes) {
	Iterator<String> it = genes.keySet().iterator();
	String[][] array_genes = new String[genes.size()][2];
	int counter = 0;

	while (it.hasNext()) {
	    String key = it.next();
	    array_genes[counter][0] = key;
	    array_genes[counter][1] = genes.get(key);
	    counter++;
	}

	return array_genes;
    }
}
