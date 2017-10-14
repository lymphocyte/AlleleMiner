package tools;

import java.io.BufferedReader;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.IOException;
import java.util.HashMap;

public class genmoeSampleFileParser {
    private String file_path = "";
    private HashMap<String, String> sample_with_origin = new HashMap<String, String>();

    public genmoeSampleFileParser() {
    }

    public genmoeSampleFileParser(String file_path) {
	this.file_path = file_path;
    }

    public void read_sample_file() throws IOException {
	BufferedReader br = new BufferedReader(new FileReader(this.file_path));
	String line = br.readLine();
	line = br.readLine();
	while (line != null) {

	    String[] tmp = line.split(",");
	    this.sample_with_origin.put(tmp[0], tmp[1]);
	    line = br.readLine();
	}
    }

    public String GetOriginByID(String id) {
	if (this.sample_with_origin.containsKey(id)) {
	    return this.sample_with_origin.get(id);
	} else {
	    System.out.println("This id doesn't exist!");
	    return "";
	}
    }

    public HashMap<String, String> getSample_map() {
	return this.sample_with_origin;
    }
}
