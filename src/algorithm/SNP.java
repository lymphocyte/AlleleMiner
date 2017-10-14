package algorithm;

public class SNP {
    private int[] location;
    private String ref_gene;
    private String alt_gene;

    public SNP(int location, String ref_gene, String alt_gene) {
	/*
	 * If there are several nucleotides on the gene changed.
	 */
	this.ref_gene = ref_gene;
	this.alt_gene = alt_gene;
	if (ref_gene.length() > 1) {
	    this.location = new int[ref_gene.length()];
	    for (int i = 0; i < this.location.length; i++) {
		this.location[i] = location;
		location++;
	    }
	}
    }
    
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
