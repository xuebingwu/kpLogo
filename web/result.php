<!DOCTYPE HTML PUBLIC "-//W3C//DTD HTML 4.01 Transitional//EN">

<html>

<head>

<meta http-equiv="Content-Type" content="text/html; charset=gb2312"><title>PKA:positional K-mer analysis</title>

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
height:30px;
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

a:link {
    text-decoration: none;
}

</style>

</head>

  <body style="margin:0 auto; width: 970px; text-align: center" >

<div style="font-family: Helvetica,Arial,sans-serif;text-align: left;" id="title"> <br>
<h1><span class="title"> PKA: Positional K-mer Analysis</span><br> </h1>
</div>

<div style=" border-top: 8px solid DimGray; width: 970px"></div>
<div style=" border-top: 8px solid White; width: 970px" ></div>


<table style="text-align: left; width: 436px; height: 32px; margin-left: auto; margin-right: 0px;" border="0" cellpadding="0" cellspacing="0">
  <tbody>
    <tr>
      <td style="vertical-align: top; width: 69px;">
      <form action="http://pka.wi.mit.edu/index.html"> <input value="Home" type="submit" style="color: GreenYellow" ></form>
      </td>
      <td style="vertical-align: top; width: 69px;">
      <form action="http://pka.wi.mit.edu/submit.html"> <input value="Server" type="submit"></form>
      </td>
      <td style="vertical-align: top; width: 0px;">
        <form action="http://pka.wi.mit.edu/manual.html#install-pka-locally"> <input value="Code" type="submit"></form>
        </td>
        <td style="vertical-align: top; width: 69px;">
        <form action="http://pka.wi.mit.edu/manual.html"> <input value="Manual" type="submit"></form>
        </td>
        <td style="vertical-align: top; width: 69px;">
        <form action="http://pka.wi.mit.edu/manual.html#examples"> <input value="Examples" type="submit"></form>
      </td>
      <td style="vertical-align: top; width: 0px;">
      <form action="http://pka.wi.mit.edu/feedback.html"> <input value="Feedback" type="submit"></form>
      </td>
    </tr>
  </tbody>
</table>

<div style="font-family: Helvetica,Arial,sans-serif;text-align: left;" > <br>


<?php 

// load job info
$handle = @fopen("jobinfo.txt", "r");
$folder = str_replace(array("\r", "\n"), '', fgets($handle, 4096));
$email = str_replace(array("\r", "\n"), '', fgets($handle, 4096));
$jobname = str_replace(array("\r", "\n"), '', fgets($handle, 4096));
$command = str_replace(array("\r", "\n"), '', fgets($handle, 4096));
$jobID = str_replace(array("\r", "\n"), '', fgets($handle, 4096));
fclose($handle);

// if exit due to error
if(file_exists("exit_with_error") == true){
	echo "<font color='red'>Exit due to error!</font><br>";
	sleep(1);
    echo nl2br(file_get_contents( "log" ));
	exit();
}

if(file_exists("./pka.output.most.significant.each.position.png") == false){
    header("refresh: 1;");
}

if (file_exists('./pka.output.most.significant.each.position.png')) {
	if(file_exists('$jobID.tar.gz') == false){
		exec("tar zcvf $jobID.tar.gz *");
	}
	echo "Your output figures and data can be accessed here:  <a href=\"./$jobID.tar.gz\"> Download $jobID.tar.gz </a> (<font color=\"red\">to be removed after 72 hours!</font>) </a> <br>";
}

	echo "<div style=\"font-family: Helvetica,Arial,sans-serif;\" id=\"figures\"> <br> ";

    echo "<table border=\"0\" cellpadding=\"0\" cellspacing=\"0\" border-spacing=\"0\" border-collapse=\"collapse\" display=\"block\" >";

if(file_exists('./pka.output.freq.png')){
    echo "<tr>";
	echo "<td><strong>Frequency Logo</strong><br> <small> <a href=\"./pka.output.freq.eps\">EPS</a> | <a href=\"./pka.output.freq.pdf\">PDF</a> | <a href=\"./pka.output.freq.png\">PNG</a></td><td><a href=\"./pka.output.freq.png\"> </small> <img src=\"./pka.output.freq.png\" style=\"display:block\" width=\"100%\"> </a></td>";
    echo "</tr>";
}
if(file_exists('./pka.output.info.png')){
    echo "<tr>";
	echo "<td><strong>Information Content Logo</strong> <br> <small><a href=\"./pka.output.info.eps\">EPS</a> | <a href=\"./pka.output.info.pdf\">PDF</a> | <a href=\"./pka.output.info.png\">PNG</a></td><td><a href=\"./pka.output.info.png\"> </small> <img src=\"./pka.output.info.png\" width=\"100%\" style=\"display:block\" > </figure> </a></td> ";
    echo "</tr>";
}
if (file_exists('./pka.output.png')) {
    echo "<tr>";
		echo "<td><strong>Probability Logo</strong><br> <small> <a href=\"./pka.output.eps\">EPS</a> | <a href=\"./pka.output.pdf\">PDF</a> | <a href=\"./pka.output.png\">PNG</a> </td><td><a href=\"./pka.output.png\"> </small> <img src=\"./pka.output.png\" width=\"100%\" style=\"display:block\" ></a></td>";
    echo "</tr>";
}
if (file_exists('./pka.output.most.significant.each.position.png')) {
    echo "<tr>";
	echo "<td><strong>K-mer Logo</strong><br> <small> <a href=\"./pka.output.most.significant.each.position.eps\">EPS</a> | <a href=\"./pka.output.most.significant.each.position.pdf\">PDF</a> | <a href=\"./pka.output.most.significant.each.position.png\">PNG</a></td><td><a href=\"./pka.output.most.significant.each.position.png\"> </small><img src=\"./pka.output.most.significant.each.position.png\" style=\"display:block\" width=\"100%\"> </a></td>";
    echo "</tr>";
}
    echo "</table>";
	echo "</div>";


if (file_exists('./pka.output.most.significant.each.position.png')) {	
	echo "<h3>Tabular output </h3>";
	echo "&nbsp;&nbsp;&nbsp;&nbsp;the most significant kmer at each position : <a href=\"./pka.output.most.significant.each.position.txt\">TXT</a> <br>";
	echo "&nbsp;&nbsp;&nbsp;&nbsp;all significant kmer : <a href=\"./pka.output.pass.p.cutoff.txt\">TXT</a> <br>";
} else {

    header("refresh: 1;");

    $url = "http://$_SERVER[HTTP_HOST]$_SERVER[REQUEST_URI]";

	echo "<h2><font color=Red>Probability logo and Kmer logo will appear here once generated...</font></h2>";

    echo "<strong>Your job $jobID ($jobname) has been submitted and results will be available here for 72 hours (this page):</strong> <br><br> <a href=\"$url\">$url</a> <br><br>";

    //$list = listdir_by_date('./');
    //foreach ($list as $file)
    //{
    //  echo "$file <br>";
    //}

    // email reminder
    if (file_exists("submission_notification_not_sent") )
    {
        $subject = "PKA job submitted: $jobID ($jobname)";
        $message = "Your job $jobID ($jobname) has been submitted and results will be available here for 72 hours: \r\n\r\n $url";
        $message = wordwrap($message, 70, "\r\n");
        mail($email, $subject, $message);
        exec("rm submission_notification_not_sent");
    } elseif ($email != "xxx@yyy.com" )
    {
        echo "An email has been sent to $email. Another email will be sent once results are available.<br><br>";
    }

}

    echo "<h3> Job ID & Name </h3>";
    echo "&nbsp;&nbsp;&nbsp;&nbsp;$jobID ($jobname)<br>";

    echo "<h3> Commandline & Screen Output </h3>";
    echo "$command<br><br>";
    echo nl2br(file_get_contents( "log" ));

    echo "<h3> Input </h3>";
    echo "&nbsp;&nbsp;&nbsp;&nbsp;pka.input.txt : <a href=\"./pka.input.txt\">txt</a> <br>";

	
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
