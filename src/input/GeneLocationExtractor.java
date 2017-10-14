package input;

import java.io.ByteArrayOutputStream;
import java.io.IOException;
import java.io.PrintStream;
import java.util.Properties;

import javax.xml.parsers.ParserConfigurationException;
import javax.xml.parsers.SAXParser;
import javax.xml.parsers.SAXParserFactory;

import org.xml.sax.SAXException;

import tools.DASHandler;

/*
 * By given the chromosome, start index and end index, GeneLocationExtractor returns the DNA sequence.
 * 
 */
public class GeneLocationExtractor {
    private String chr;
    private int start_index;
    private int end_index;
    private String result_seq;
    private String DNA_database_prefix;

    GeneLocationExtractor(String chr, int start_index, int end_index, Properties props)
	    throws ParserConfigurationException, SAXException, IOException {
	this.chr = chr;
	this.start_index = start_index;
	this.end_index = end_index;
	this.DNA_database_prefix = props.getProperty("DNA_database_prefix");
	//GetDNAseq();
    }

    public String GetDNAseq() throws ParserConfigurationException, SAXException{
	// Create a stream to hold the output
	DASHandler dh = new DASHandler();
	ByteArrayOutputStream baos = new ByteArrayOutputStream();
	PrintStream ps = new PrintStream(baos);
	// IMPORTANT: Save the old System.out!
	PrintStream old = System.out;
	// Tell Java to use your special stream
	System.setOut(ps);
	SAXParserFactory f = SAXParserFactory.newInstance();
	SAXParser parser = f.newSAXParser();
	//System.out.println(DNA_database_prefix + this.chr + ":" + this.start_index + "," + this.end_index);
	try {
	    parser.parse(DNA_database_prefix + this.chr + ":" + this.start_index + "," + this.end_index, dh);
	    System.out.flush();
	} catch (IOException | SAXException e) {
	    // TODO Auto-generated catch block
	    System.setOut(old);
	    System.out.println(DNA_database_prefix + this.chr + ":" + this.start_index + "," + this.end_index);
	    System.out.println("This gene location doesn't exist!");
	    System.out.println(baos.toString());
	    System.exit(0);
	    
	}

	// Put things back
	
	System.setOut(old);
	// Show what happened
	result_seq = dh.GetDNAseq().toUpperCase();
	//System.out.println(result_seq);
	return result_seq;
    }
}
