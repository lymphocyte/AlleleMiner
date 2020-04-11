package algorithm; //create this file in algorithm package

import java.util.ArrayList; //include a library/built-in class named java.util.ArrayList

public class Haplotype {//declare a public class Haplotype

    //private data memebers
    private int start_index;
    private int end_index;
    private String ref_seq;
    private String allele;
    private String genename;

    //3 array list private data members and defualt size is 0 of all
    private ArrayList<String> ref_nt_factory = new ArrayList<String>();
    private ArrayList<String> alt_nt_factory = new ArrayList<String>();
    private ArrayList<String> locations = new ArrayList<String>();// Locations
								  // are not
								  // necessarily
								  // just one
								  // index, can
								  // be from
								  // index1-index2;
    //private data members
    private String ID;
    private String Sample_Origin;
    private String hap_seq;

    //public parameterized constructor
    public Haplotype(String ref_seq, int start_index, int end_index, String genename) {
	this.ref_seq = ref_seq; //this is used to refer the data member of the class
	this.start_index = start_index;
	this.end_index = end_index;
	this.genename = genename;
    }

    //2 setter getter for private data members ID and Origin
    public void setID(String Id) {
	this.ID = Id;
    }

    public String getID() {
	return this.ID;
    }

    public void setOrigin(String origin) {
	this.Sample_Origin = origin;
    }

    public String getOrigin() {
	return this.Sample_Origin;
    }

    //function used to set/add values in the 3 array lists defined as data members
    public void AddSNP(String ref_seq, String alt_seq, int location) {
	ref_nt_factory.add(ref_seq);
	alt_nt_factory.add(alt_seq);
	// if the ref nt is not only one, so they have several locations.

	if (ref_seq.length() > 1) {
	    locations.add(location + "-" + (location + ref_seq.length()));
	} else {
	    locations.add(location + "");
	}
    }

    //setter for ID
    public void SetID(String id) {
	this.ID = id;
    }

    public void GenerateNewHaplotype() {

	StringBuilder new_allele_builder = new StringBuilder(this.ref_seq);  //an object of StringBuilder class just like String class 

	for (int i = 0; i < locations.size(); i++) { //retrieve the size of locations array list
	    String tmp_loc = locations.get(i); //get the locations value from specific index i
	    if (tmp_loc.contains("-")) { //if tmp_loc contains "-" then return true
		if (Integer.parseInt(tmp_loc.split("-")[0]) > new_allele_builder.length()) {
		/*
		-> parseInt: it will convert the string into integer
		->split: it will split the string on the basis of any symbol
		->tmp_loc.split("-")[0]: it will split the tmp_loc on the basis of "-" and retrieve only index 0 value
		->.length(): it will return the size of StringBuilder object new_allele_builder
		*/
		    this.hap_seq = new_allele_builder.toString(); //toString(): it will convert it in string and return
		    return;
		}
		int start_index = Integer.parseInt(tmp_loc.split("-")[0]);
		int end_index = Integer.parseInt(tmp_loc.split("-")[1]);
		new_allele_builder.replace(start_index, end_index, alt_nt_factory.get(i)); //it will replace contents of new_allele_builder from start_index to end_index by alt_nt_factory.get(i)
	    } else {
		// System.out.println(new_allele_builder.charAt(Integer.parseInt(tmp_loc.split("-")[0])));
		// System.out.println(Integer.parseInt(tmp_loc.split("-")[0])+"
		// "+new_allele_builder.length());
		if (Integer.parseInt(tmp_loc.split("-")[0]) > new_allele_builder.length()) {
		    this.hap_seq = new_allele_builder.toString();
		    return;
		}
		new_allele_builder.replace(Integer.parseInt(tmp_loc.split("-")[0]),
			Integer.parseInt(tmp_loc.split("-")[0]) + 1, alt_nt_factory.get(i));
	    }
	}

	this.hap_seq = new_allele_builder.toString();
	this.allele = this.hap_seq;

    }
    //again setter getters
    public String getHap_seq() {
	return this.hap_seq;
    }

    public String getAllele() {
	return this.allele;
    }

    public void SetRefSeq(String ref_seq) {
	this.ref_seq = ref_seq;
    }

    public String getGeneName() {
	return this.genename;
    }

}
