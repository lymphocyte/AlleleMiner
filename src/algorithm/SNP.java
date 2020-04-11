package algorithm; //file is created in a package named algorithm

public class SNP {  { //class name
    //3 data members
    private int[] location;
    private String ref_gene;
    private String alt_gene;

    public SNP(int location, String ref_gene, String alt_gene) {
	/*
	 * If there are several nucleotides on the gene changed.
	 */
	this.ref_gene = ref_gene; //this is used to refer the attribute of the class i.e. SNP
	this.alt_gene = alt_gene;
	if (ref_gene.length() > 1) { //.length() is the function of String to determine its length
	    this.location = new int[ref_gene.length()];
	    for (int i = 0; i < this.location.length; i++) {
		this.location[i] = location;
		location++;
	    }
	}
    }
    //3 getter functions of private data members
    public int[] getLocation()
    {
	return this.location;
    }
    
    public String get_ref_gene()
    {
	return this.ref_gene;
    }
    
    public String get_alt_gene()
    {
	return this.alt_gene;
    }
}
