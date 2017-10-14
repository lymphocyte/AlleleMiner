package tools;

import javax.xml.parsers.SAXParser;
import javax.xml.parsers.SAXParserFactory;

import org.xml.sax.Attributes;
import org.xml.sax.SAXException;
import org.xml.sax.helpers.DefaultHandler;


public class DASHandler extends DefaultHandler {
    private boolean inDNA = false;
    private String result = "";
    @Override
    public void startElement(String uri, String localName, String qName, Attributes attributes) throws SAXException {
	inDNA = (qName.equals("DNA"));
    }

    @Override
    public void endElement(String uri, String localName, String qName) throws SAXException {
	inDNA = false;
    }

    @Override
    public void characters(char[] ch, int start, int length) throws SAXException {
	
	if (inDNA)
	{
	    result +=new String(ch, start, length).replace("\n", "");
	}
    }
    
    public String GetDNAseq()
    {
	return this.result;
    }
    
}