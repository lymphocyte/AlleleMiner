package tools;

public class ReverseCompliment {

	private String seq;
	
	public  ReverseCompliment(String seq)
	{
		this.seq = seq;
	}
	
	public String getReverse()
	{
		String seq = new StringBuilder(this.seq).reverse().toString();
		return seq;
	}
	
	public String getReverseCompliment()
	{
		String reversed_seq = getReverse();
		StringBuilder tmpStr = new StringBuilder();
		for(int i=0;i<reversed_seq.length();i++)
		{
			if(reversed_seq.charAt(i)=="A".charAt(0))
			{
				tmpStr.append("T");
			}else if(reversed_seq.charAt(i)=="T".charAt(0))
			{
				tmpStr.append("A");
			}else if(reversed_seq.charAt(i)=="C".charAt(0))
			{
				tmpStr.append("G");
			}else if(reversed_seq.charAt(i)=="G".charAt(0))
			{
				tmpStr.append("C");
			}else
			{
				tmpStr.append(reversed_seq.charAt(i));
			}
		}
		
		return tmpStr.toString();
	}
}
