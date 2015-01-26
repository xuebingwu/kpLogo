<?php
// Start the buffering //
ob_start();
?>

<!DOCTYPE HTML PUBLIC "-//W3C//DTD HTML 4.01 Transitional//EN">
<html><head>


  
  <meta http-equiv="Content-Type" content="text/html; charset=gb2312"><title>PKA:positional kmer analysis</title>
  

  
  
  <style>
input[type="submit"]{
/* change these properties to whatever you want */
background-color: #555;
color: #fff;
border-radius: 10px;
font-size:1em;
}
  </style>
  
  <style>
input[type="file"]{
/* change these properties to whatever you want */
background-color: #555;
color: #fff;
border-radius: 10px;
font-size:1em;
}
  </style>
  
  <style>
input[type="text"]{
/* change these properties to whatever you want */
background-color: #aaa;
color: #fff;
}
  </style></head><body style="margin-left: 101px; width: 970px;">
<div style="font-family: Helvetica,Arial,sans-serif;" id="title"> <br>
<h1><span class="title">
PKA: Positional Kmer Analysis</span><br>
</h1>
<span class="title"></span>
<hr style="width: 100%; height: 8px;" noshade="noshade"></div>

</body></html>

<?php // get values of all variables

$foreground=$_POST['foreground'];	// foreground sequences

$forgroundfile=$_POST['forgroundfile'];	// forgroundfile

$alphabet=$_POST['alphabet'];		// alphabet: DNA, protein, other

$kmer_length = $_POST['kmer_length'];

$kmer_upto = $_POST['kmer_upto'];

$min_shift = $_POST['min_shift'];

$max_shift = $_POST['max_shift'];

$degenerate = $_POST['degenerate'];

$degenerate_alphabet = $_POST['degenerate_alphabet'];

$background = $_POST['background']; // markov_foreground, shuffle, backgroundfromfile, markov_background

$markov_foreground_order = $_POST['markov_foreground_order'];

$shuffle_n = $_POST['shuffle_n'];

$shuffle_m = $_POST['shuffle_m'];

$backgroundfile = $_POST['backgroundfile'];

$output = $_POST['output'];


// make a folder with randome name for each job
$randomString = substr(str_shuffle(md5(time())),0,20);
 
$tmpfolder = "./files_to_be_removed_72_hrs_after_creation/".$randomString;

$oldmask = umask(0);
if (!mkdir($tmpfolder, 0777, true)) {
    die('Failed to create folders...');
}

// enter this folder
chdir($tmpfolder);

// input
if (strlen($foreground) > 0) // sequence pasted, save to file
{
	$h = fopen("PKA.input.txt", 'w');
	fwrite($h, $foreground);
	fclose($h); 
}
else
{
	if(move_uploaded_file($_FILES['forgroundfile']['tmp_name'], "PKA.input.txt")) 
	{
		//echo "The file ".basename( $_FILES['uploadedfile']['name'])." has been uploaded";
	}
	else
	{
		echo "<font color='red'>Please paste or upload your input!</font>";
		exit();
	}
}


$k = "-k";
if ($kmer_upto == "on")
{
	$k = "-upto";
}

$command = "../../PKA PKA.input.txt -o PKA";




$result = exec($command);


//$myFile = "result.html"; // or .php   
//$fh = fopen($myFile, 'w'); // or die("error");  
//$stringData = "you html code php code goes here";   
//fwrite($fh, $stringData);
//fclose($fh);

umask($oldmask);

$path_parts = pathinfo($_SERVER['REQUEST_URI']);
$path = $path_parts['dirname'];

$url_folder = "http" . ($_SERVER['HTTPS'] ? 's' : '') . "://{$_SERVER['HTTP_HOST']}" . $path . "/files_to_be_removed_72_hrs_after_creation/".$randomString;

echo "Your output figures and data can be accessed here:  <a href=\"$url_folder/all.tar.gz\"> Download </a> (<font color=\"red\">to be removed after 72 hours!</font>) <br> This page: <a href=\"$url_folder/result.html\"> $url_folder/output.html </a> <br>";

echo "<div style=\"font-family: Helvetica,Arial,sans-serif;\" id=\"figures\"> <br> ";
	
echo "<a href=\"$url_folder/PKA.freq.png\"> <figure> <img src=\"$url_folder/PKA.freq.png\" width=\"400\"> </figure> </a>";
echo "<a href=\"$url_folder/PKA.info.png\"> <figure> <img src=\"$url_folder/PKA.info.png\" width=\"400\"> </figure> </a> ";
echo "<a href=\"$url_folder/PKA.png\"> <figure> <img src=\"$url_folder/PKA.png\" width=\"400\"> </figure> </a>";
echo "<a href=\"$url_folder/PKA.most.significant.each.position.png\"> <figure> <img src=\"$url_folder/PKA.most.significant.each.position.png\" width=\"400\"> </figure> </a>";

echo "</div>";

echo "<div style=\"font-family: Helvetica,Arial,sans-serif;\" id=\"info\"> <br> ";

echo "<strong> commandline </strong> <br>";

echo $command."<br>";

echo "</div>";


echo "Input: <a href=\"$url_folder/PKA.input.txt\">download</a> <br>";

echo "Output:<br>";
echo "\t<a href=\"$url_folder/PKA.freq.ps\">frequency logo</a> <br>";
echo "\t<a href=\"$url_folder/PKA.info.ps\">information content logo</a> <br>";
echo "\t<a href=\"$url_folder/PKA.ps\">significance logo</a> <br>";
echo "\t<a href=\"$url_folder/PKA.most.significant.each.position.ps\">kmer logo</a> <br>";
echo "\t<a href=\"$url_folder/PKA.most.significant.each.position.txt\">most.significant.kmer.at.each.position</a> <br>";
echo "\t<a href=\"$url_folder/PKA.pass.p.cutoff.txt\">output data file</a> <br>";

?>

<?php
// Get the content that is in the buffer and put it in your file //
file_put_contents('result.html', ob_get_contents());
?>