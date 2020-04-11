package input; //package input

//include java packages
import java.io.ByteArrayOutputStream;
import java.io.IOException;
import java.io.PrintStream;
import java.util.Properties;

//include javax packages
import javax.xml.parsers.ParserConfigurationException;
import javax.xml.parsers.SAXParser;
import javax.xml.parsers.SAXParserFactory;

//include SAXException
import org.xml.sax.SAXException; //used for SAX exception, SAX mean Simple Api for XML

//including some other package
import tools.DASHandler;

/*
 * By given the chromosome, start index and end index, GeneLocationExtractor returns the DNA sequence.
 * 
 */
public class GeneLocationExtractor { //main class
    private String chr;
    private int start_index;
    private int end_index;
    private String result_seq;
    private String DNA_database_prefix;

    //parameterized constructor
    GeneLocationExtractor(String chr, int start_index, int end_index, Properties props)
	    throws ParserConfigurationException, SAXException, IOException {
	this.chr = chr;
	this.start_index = start_index;
	this.end_index = end_index;
	this.DNA_database_prefix = props.getProperty("DNA_database_prefix");//getProperty: search the value against the given key in paramter
	//GetDNAseq();
    }

    public String GetDNAseq() throws ParserConfigurationException, SAXException{
	// Create a stream to hold the output
	DASHandler dh = new DASHandler();
	ByteArrayOutputStream baos = new ByteArrayOutputStream();//This class implements an output stream in which the data is written into a byte array. The buffer automatically grows as data is written to it.
	PrintStream ps = new PrintStream(baos); //printStream for writing/printing
	// IMPORTANT: Save the old System.out!
	PrintStream old = System.out; //Stream for writing to console - standard output stream
	// Tell Java to use your special stream
	System.setOut(ps); //Now, from console it is changed to printStream object ps, which is basically our bytearray
	SAXParserFactory f = SAXParserFactory.newInstance(); //Defines a factory API that enables applications to configure and obtain a SAX based parser to parse XML documents.
	SAXParser parser = f.newSAXParser(); //Creates a new instance of a SAXParser using the currently configured factory parameters.
	//System.out.println(DNA_database_prefix + this.chr + ":" + this.start_index + "," + this.end_index);
	try {
	    parser.parse(DNA_database_prefix + this.chr + ":" + this.start_index + "," + this.end_index, dh); //see at end
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
/*
-> public void parse(String uri, DefaultHandler dh) throws SAXException, IOException

Parse the content described by the giving Uniform Resource Identifier (URI) as XML using the specified DefaultHandler.

Parameters:
uri - The location of the content to be parsed.
dh - The SAX DefaultHandler to use.

Throws:
IllegalArgumentException - If the uri is null.
IOException - If any IO errors occur.
SAXException - If any SAX errors occur during processing.
*/
