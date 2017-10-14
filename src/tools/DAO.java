package tools;

import java.sql.Connection;
import java.sql.DriverManager;
import java.sql.ResultSet;
import java.sql.SQLException;
import java.sql.Statement;

import javax.sql.PooledConnection;

/*
 *  Database Access Object
 */
public class DAO {
    
    String username;
    String host;
    String pwd;
    //PooledConnection pool;
    Connection con;

    
    
    
    public DAO(String host, String username, String pwd) throws SQLException
    {
	this.username = username;
	this.host = host;
	this.pwd = pwd;
	InianizeConnection();
    }
    
    private void InianizeConnection() throws SQLException
    {	  
	//Connection con2 =  DriverManager.getConnection(host,username,pwd);
         con = DriverManager.getConnection("jdbc:mysql://genome-mysql.cse.ucsc.edu/hg19","genome", "");

	//con = con2;
    }
    
    public ResultSet search(String query) throws SQLException
    {
	 ResultSet result;
	 String sql_seq = query;
	 Statement pstmt= con.createStatement();
	 try {
	     result = pstmt.executeQuery(sql_seq);
	} catch (SQLException e) {
	    // TODO Auto-generated catch block
	    //e.printStackTrace();
	    return null;
	}
	 
	return result;
    }
    
    public void insert()
    {
	
    }
    
    public void delete()
    {
	
    }
    
    public void modify()
    {
	
    }
}
