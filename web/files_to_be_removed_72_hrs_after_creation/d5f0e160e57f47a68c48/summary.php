<!DOCTYPE HTML PUBLIC "-//W3C//DTD HTML 4.01 Transitional//EN">

<html>

<head>

<meta http-equiv="Content-Type" content="text/html; charset=gb2312"><title>PKA:positional kmer analysis</title>

<style>

input {
text-align: center;
}

input[type="submit"]{
/* change these properties to whatever you want */
background-color: #555;
color: #fff;
border-radius: 10px;
font-size:1em;
}

input[type="file"]{
/* change these properties to whatever you want */
background-color: #aaa;
color: #fff;
border-radius: 10px;
/*font-size:1em;*/
}
  
input[type="text"]{
/* change these properties to whatever you want */
background-color: #aaa;
color: #fff;
}

</style>

</head>

  <body style="margin:0 auto; width: 970px; text-align: center" >
<div style="font-family: Helvetica,Arial,sans-serif;text-align: left;" id="title"> <br>
	
<h1><span class="title"> PKA: Positional Kmer Analysis</span><br> </h1>
<hr style="width: 100%; height: 8px;" noshade="noshade"></div>

<table style="text-align: left; width: 436px; height: 32px; margin-left: auto; margin-right: 0px;" border="0" cellpadding="0" cellspacing="0">
  <tbody>
    <tr>
      <td style="vertical-align: top; width: 69px;">
      <form action="../../index.html"> <input value="Home" type="submit"></form>
      </td>
      <td style="vertical-align: top; width: 69px;">
      <form action="../../submit.html"> <input value="Server" type="submit"></form>
      </td>
      <td style="vertical-align: top; width: 0px;">
      <form action="../../code.html"> <input value="Code" type="submit"></form>
      </td>
      <td style="vertical-align: top; width: 69px;">
      <form action="../../manual.html"> <input value="Manual" type="submit"></form>
      </td>
      <td style="vertical-align: top; width: 69px;">
      <form action="../../examples.html"> <input value="Examples" type="submit"></form>
      </td>
	  <td style="vertical-align: top; width: 0px;">
	  <form action="../../feedback.html"> <input value="Feedback" type="submit"></form>
	  </td>
    </tr>
  </tbody>
</table>

<div style="font-family: Helvetica,Arial,sans-serif;text-align: left;" > <br>


<?php 

session_start();




if (file_exists('./pka.output.most.significant.each.position.png')) {

	exec("tar zcvf pka.output.result.tar.gz *");

	echo "Your output figures and data can be accessed here:  <a href=\"./pka.output.result.tar.gz\"> Download </a> (<font color=\"red\">to be removed after 72 hours!</font>) </a> <br>";

	echo "<div style=\"font-family: Helvetica,Arial,sans-serif;\" id=\"figures\"> <br> ";
	
	echo "<a href=\"./pka.output.freq.png\"> <figure> <img src=\"./pka.output.freq.png\" width=\"400\"> </figure> </a>";
	echo "<a href=\"./pka.output.info.png\"> <figure> <img src=\"./pka.output.info.png\" width=\"400\"> </figure> </a> ";
	if (file_exists('./pka.output.png')) {
		echo "<a href=\"./pka.output.png\"> <figure> <img src=\"./pka.output.png\" width=\"400\"> </figure> </a>";
	}
	echo "<a href=\"./pka.output.most.significant.each.position.png\"> <figure> <img src=\"./pka.output.most.significant.each.position.png\" width=\"400\"> </figure> </a>";

	echo "</div>";


    echo "<h3> Job name </h3>";
    echo "&nbsp;&nbsp;&nbsp;&nbsp;" . $_SESSION['jobname'] . "<br>";

	echo "<h3> Commandline </h3>";
	echo "&nbsp;&nbsp;&nbsp;&nbsp;" . $_SESSION['command'] . "<br>";

	echo "<h3> Input </h3>";
	echo "&nbsp;&nbsp;&nbsp;&nbsp;pka.input.txt : <a href=\"./pka.input.txt\">txt</a> <br>";
	
	echo "<h3>Logo output </h3>";
	echo "&nbsp;&nbsp;&nbsp;&nbsp;frequency logo : <a href=\"./pka.output.freq.ps\">ps</a> <a href=\"./pka.output.freq.pdf\">pdf</a> <a href=\"./pka.output.freq.png\">png</a><br>";
	echo "&nbsp;&nbsp;&nbsp;&nbsp;information content logo : <a href=\"./pka.output.info.ps\">ps</a> <a href=\"./pka.output.info.pdf\">pdf</a> <a href=\"./pka.output.info.png\">png</a><br>";
	echo "&nbsp;&nbsp;&nbsp;&nbsp;probability logo : <a href=\"./pka.output.ps\">ps</a> <a href=\"./pka.output.pdf\">pdf</a> <a href=\"./pka.output.png\">png</a> <br>";
	echo "&nbsp;&nbsp;&nbsp;&nbsp;kmer logo : <a href=\"./pka.output.most.significant.each.position.ps\">ps</a> <a href=\"./pka.output.most.significant.each.position.pdf\">pdf</a> <a href=\"./pka.output.most.significant.each.position.png\">png</a> <br>";
	echo "<br>";
	echo "<h3>Tabular output </h3>";
	echo "&nbsp;&nbsp;&nbsp;&nbsp;the most significant kmer at each position : <a href=\"./pka.output.most.significant.each.position.txt\">txt</a> <br>";
	echo "&nbsp;&nbsp;&nbsp;&nbsp;raw data : <a href=\"./pka.output.pass.p.cutoff.txt\">txt</a> <br>";

	// email
	$email = $_SESSION['email'];
	if ($email != "xxx@yyy.com")
	{
		//echo "$email <br>";
		$url = "http://$_SERVER[HTTP_HOST]$_SERVER[REQUEST_URI]";
		//echo "$url <br>";
		$message = "Your PKA job '" . $_SESSION['jobname'] . "' has been finished and results are available here within 72 hours: \r\n\r\n $url ";
		// In case any of our lines are larger than 70 characters, we should use wordwrap()
		$message = wordwrap($message, 70, "\r\n");
        $subject = 'PKA results available: ' . $_SESSION['jobname'];
		mail($email, $subject, $message);
		$_SESSION['email'] = "xxx@yyy.com"; // avoid send emails every time this page is loaded, only send email once
		
	}

}
else{
	echo " <br>Your job has been submitted and results will be available on this page. Typically it takes a few seconds. <br><br><br>";
	
	//$list = listdir_by_date('./');
	//foreach ($list as $file)
	//{
	//	echo "$file <br>";
	//}

	header("refresh: 1;");

}


	
	function listdir_by_date($path){
	$dir = opendir($path);
	$list = array();
	while($file = readdir($dir)){
	if ($file != '.' and $file != '..'){
	// add the filename, to be sure not to
	// overwrite a array key
	$ctime = filectime($path . $file) . ',' . $file;
	$list[$ctime] = $file;
	}
	}
	closedir($dir);
	ksort($list); // krsort
	return $list;
	}
	

?>

	   <hr> <p align="center">Copyright (c) 2014 Xuebing Wu <br><br><br>
		   
	   </div>

</body></html>
