<!DOCTYPE HTML PUBLIC "-//W3C//DTD HTML 4.01 Transitional//EN">

<html>

<head>

<meta http-equiv="Content-Type" content="text/html; charset=gb2312"><title>kpLogo: k-mer probability logo</title>

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
<h1><span class="title"> <i>k</i>pLogo:  <i>k</i>-mer probability logo</span><br> </h1>
</div>

<div style=" border-top: 8px solid DimGray; width: 970px"></div>
<div style=" border-top: 8px solid White; width: 970px" ></div>

<table style="text-align: left; width: 436px; height: 32px; margin-left: auto; margin-right: 0px;" border="0" cellpadding="0" cellspacing="0">
  <tbody>
    <tr>
      <td style="vertical-align: top; width: 69px;">
      <form action="http://kplogo.wi.mit.edu/index.html"> <input value="Home" type="submit" style="color: GreenYellow" ></form>
      </td>
      <td style="vertical-align: top; width: 69px;">
      <form action="http://kplogo.wi.mit.edu/submit.html"> <input value="Server" type="submit"></form>
      </td>
      <td style="vertical-align: top; width: 0px;">
        <form action="http://kplogo.wi.mit.edu/manual.html#install-kpLogo-locally"> <input value="Code" type="submit"></form>
        </td>
        <td style="vertical-align: top; width: 69px;">
        <form action="http://kplogo.wi.mit.edu/manual.html"> <input value="Manual" type="submit"></form>
        </td>
        <td style="vertical-align: top; width: 69px;">
        <form action="http://kplogo.wi.mit.edu/manual.html#examples"> <input value="Examples" type="submit"></form>
      </td>
	  <td style="vertical-align: top; width: 0px;">
	  <form action="http://kplogo.wi.mit.edu/feedback.html"> <input value="Feedback" type="submit"></form>
	  </td>
    </tr>
  </tbody>
</table>

<div style="font-family: Helvetica,Arial,sans-serif;text-align: left;" > <br>


<?php 

// load job info
$handle = @fopen("job.info", "r");
$folder = str_replace(array("\r", "\n"), '', fgets($handle, 4096));
$jobname = str_replace(array("\r", "\n"), '', fgets($handle, 4096));
$command = str_replace(array("\r", "\n"), '', fgets($handle, 4096));
$jobID = str_replace(array("\r", "\n"), '', fgets($handle, 4096));
fclose($handle);

$email = "";

// if exit due to error
if(file_exists("exit_with_error") == true){
	echo "<font color='red'>Exit due to error!</font><br>";
	sleep(1);
    echo nl2br(file_get_contents( "log" ));
	exec("rm *");
	exit();
}

if(file_exists("./done") == false){
    header("refresh: 1;");
} else {
	if(file_exists('$jobID.tar.gz') == false){
		exec("tar zcvf $jobID.tar.gz *");
	}
	echo "Your output figures and data can be accessed here:  <a href=\"./$jobID.tar.gz\"> Download $jobID.tar.gz </a> (<font color=\"red\">to be removed after 10 days!</font>) </a> <br>";
	
	
	
}

	echo
	"<form action='' method='post'>
	<input type='submit' name='delete_data' value='Delete all data from the server' style='color: Yellow'  />
	</form>";

	if(isset($_POST['delete_data']))
	{
	    exec("rm ./*");
		//echo '<script type="text/javascript">alert("All data have been deleted from the server. Thanks for using kpLogo!");</script>';
		header("Location: ../../index.html");
		exit();
	}

	echo "<div style=\"font-family: Helvetica,Arial,sans-serif;\" id=\"figures\"> <br> ";

    echo "<table border=\"0\" cellpadding=\"0\" cellspacing=\"0\" border-spacing=\"0\" border-collapse=\"collapse\" display=\"block\" >";

if(file_exists('./kpLogo.output.freq.png')){
    echo "<tr>";
	echo "<td><strong>Frequency Logo</strong><br> <small> <a href=\"./kpLogo.output.freq.eps\">EPS</a> | <a href=\"./kpLogo.output.freq.pdf\">PDF</a> | <a href=\"./kpLogo.output.freq.png\">PNG</a></td><td><a href=\"./kpLogo.output.freq.png\"> </small> <img src=\"./kpLogo.output.freq.png\" style=\"display:block\" width=\"100%\"> </a></td>";
    echo "</tr>";
}
if(file_exists('./kpLogo.output.info.png')){
    echo "<tr>";
	echo "<td><strong>Information Content Logo</strong> <br> <small><a href=\"./kpLogo.output.info.eps\">EPS</a> | <a href=\"./kpLogo.output.info.pdf\">PDF</a> | <a href=\"./kpLogo.output.info.png\">PNG</a></td><td><a href=\"./kpLogo.output.info.png\"> </small> <img src=\"./kpLogo.output.info.png\" width=\"100%\" style=\"display:block\" > </figure> </a></td> ";
    echo "</tr>";
}
if (file_exists('./kpLogo.output.most.significant.each.position.png')) {
    echo "<tr>";
	echo "<td><strong><i>k</i>-mer Logo</strong><br> <small> <a href=\"./kpLogo.output.most.significant.each.position.eps\">EPS</a> | <a href=\"./kpLogo.output.most.significant.each.position.pdf\">PDF</a> | <a href=\"./kpLogo.output.most.significant.each.position.png\">PNG</a></td><td><a href=\"./kpLogo.output.most.significant.each.position.png\"> </small><img src=\"./kpLogo.output.most.significant.each.position.png\" style=\"display:block\" width=\"100%\"> </a></td>";
    echo "</tr>";
}
if (file_exists('./kpLogo.output.png')) {
    echo "<tr>";
		echo "<td><strong>Probability Logo</strong><br> <small> <a href=\"./kpLogo.output.eps\">EPS</a> | <a href=\"./kpLogo.output.pdf\">PDF</a> | <a href=\"./kpLogo.output.png\">PNG</a> </td><td><a href=\"./kpLogo.output.png\"> </small> <img src=\"./kpLogo.output.png\" width=\"100%\" style=\"display:block\" ></a></td>";
    echo "</tr>";
}
    echo "</table>";
	echo "</div>";


if (file_exists('./done')) {	
	echo "<h3>Tabular output </h3>";
	echo "&nbsp;&nbsp;&nbsp;&nbsp;the most significant kmer at each position : <a href=\"./kpLogo.output.most.significant.each.position.txt\">TXT</a> <br>";
	echo "&nbsp;&nbsp;&nbsp;&nbsp;all significant kmer : <a href=\"./kpLogo.output.pass.p.cutoff.txt\">TXT</a> <br>";
	
	if (file_exists("./result_email_not_sent")){

		$handle1 = @fopen("result_email_not_sent", "r");
		$email = str_replace(array("\r", "\n"), '', fgets($handle1, 4096));
		fclose($handle1);
		
    $subject = "kpLogo results available: $jobID ($jobname)";
    $message = "Your job $jobID ($jobname) is finished. The logo plots depicting k-mer probability, information content, frequency, and monomer probability (if k=1 included) are attached. More detailed resutls can be found here (to be removed in *** 10 days ***): \r\n\r\n http://$_SERVER[HTTP_HOST]$_SERVER[REQUEST_URI]";
    $message = wordwrap($message, 70, "\r\n");
	
	$path = "./";
	$filename = "kpLogo.output.all.pdf";
	$file = $path.$filename;
	$content = file_get_contents($file);
	$content = chunk_split(base64_encode($content));
	$uid = md5(uniqid(time()));
	$name = basename($file);
	$from_name = "kpLogo";
	$from_mail = "www-data@wi.mit.edu";
	$replyto = "wuxb07@gmail.com";

	// header
	$header = "From: ".$from_name." <".$from_mail.">\r\n";
	$header .= "Reply-To: ".$replyto."\r\n";
	$header .= "MIME-Version: 1.0\r\n";
	$header .= "Content-Type: multipart/mixed; boundary=\"".$uid."\"\r\n\r\n";

	// message & attachment
	$nmessage = "--".$uid."\r\n";
	$nmessage .= "Content-type:text/plain; charset=iso-8859-1\r\n";
	$nmessage .= "Content-Transfer-Encoding: 7bit\r\n\r\n";
	$nmessage .= $message."\r\n\r\n";
	$nmessage .= "--".$uid."\r\n";
	$nmessage .= "Content-Type: application/octet-stream; name=\"".$filename."\"\r\n";
	$nmessage .= "Content-Transfer-Encoding: base64\r\n";
	$nmessage .= "Content-Disposition: attachment; filename=\"".$filename."\"\r\n\r\n";
	$nmessage .= $content."\r\n\r\n";
	$nmessage .= "--".$uid."--";

	if (mail($email, $subject, $nmessage, $header)) {
		exec("rm ./result_email_not_sent");
	}
	
	}
} else {

    header("refresh: 1;");

    $url = "http://$_SERVER[HTTP_HOST]$_SERVER[REQUEST_URI]";

	echo "<h2><font color=Red>Probability logo and <i>k</i>-mer logo will appear here once generated...Do NOT close this page to receive an email notification</font></h2>";

    echo "<strong>Your job $jobID ($jobname) has been submitted and results will be available here for 10 days (this page):</strong> <br><br> <a href=\"$url\">$url</a> <br> You need to leave this page open to receive an email notification when the result is available.<br>";

    //$list = listdir_by_date('./');
    //foreach ($list as $file)
    //{
    //  echo "$file <br>";
    //}

    // email reminder
    if (file_exists("submission_notification_not_sent") )
    {
		$handle1 = @fopen("submission_notification_not_sent", "r");
		$email = str_replace(array("\r", "\n"), '', fgets($handle1, 4096));
		fclose($handle1);
		
        $subject = "kpLogo job submitted: $jobID ($jobname)";
        $message = "Your job $jobID ($jobname) has been submitted. The results will be available at the link below. Once available, the data will be retained on the server for 10 days, or you can delete the data from the server immediately. To receive another email when the results are available, please click the following link and keep it open: \r\n\r\n $url";
        $message = wordwrap($message, 70, "\r\n");
		if(mail($email, $subject, $message)){
	        exec("mv submission_notification_not_sent result_email_not_sent");
		} else {
			echo "<font color='red'>Failed to send email to $email! Make sure you entered a correct email address</font><br>";
		}
    }
}

    echo "<h3> Job ID & Name </h3>";
    echo "&nbsp;&nbsp;&nbsp;&nbsp;$jobID ($jobname)<br>";

    echo "<h3> Commandline & Screen Output </h3>";
    echo "$command<br><br>";
    echo nl2br(file_get_contents( "log" ));

    echo "<h3> Input </h3>";
    echo "&nbsp;&nbsp;&nbsp;&nbsp;kpLogo.input.txt : <a href=\"./kpLogo.input.txt\">txt</a> <br>";

	
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
