

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

<div style="font-family: Helvetica,Arial,sans-serif;text-align: left;" id="title"> <br>
<h1><span class="title"> PKA: Positional Kmer Analysis</span><br> </h1>
</div>

<div style=" border-top: 8px solid DimGray; width: 970px"></div>
<div style=" border-top: 8px solid White; width: 970px" ></div>


</body></html>

<?php // get values of all variables

session_start();

// $command = "../../PKA PKA.input.txt -o PKA $inputtype -seq $col_seq -weight $col_weight -alphabet $alphabet  $kmer_length -shift $shift $background -startPos $startPos";

$email=$_POST['email'];
$jobname=$_POST['jobname'];

$foregroundpaste=$_POST['foregroundpaste'];	// foreground sequences
$foregroundfile=$_POST['foregroundfile'];	// forgroundfile
$backgroundpaste=$_POST['backgroundpaste'];	// foreground sequences
$backgroundfile=$_POST['backgroundfile'];	// forgroundfile

$inputtype=$_POST['inputtype']; // (none), -ranked, -weighted

$col_seq=$_POST['col_seq'];
$col_weight=$_POST['col_weight'];



$alphabet=$_POST['alphabet'];		// alphabet: DNA, protein, other
$other_alphabet=$_POST['other_alphabet'];		// alphabet: DNA, protein, other
if($alphabet == "other")
{
	$alphabet = $other_alphabet;
}

$kmer_length = $_POST['kmer_length'];

$shift = $_POST['shift'];

$degenerate = $_POST['degenerate'];
$degenerate_alphabet = $_POST['degenerate_alphabet'];

if($degenerate == "other")
{
	if($degenerate_alphabet == "Other: enter all allowed symbols here")
	{
		echo "<font color='red'> Please enter degenerate symbols !</font>";
		exit();	
	}
	$degenerate = " -degenerate $degenerate_alphabet ";
}


$background = $_POST['background']; // 
$markov_foreground_order = $_POST['markov_foreground_order'];
$shuffle_n = $_POST['shuffle_n'];
$shuffle_m = $_POST['shuffle_m'];

if($background == "markov_foreground"){
	$background = " -markov $markov_foreground_order ";
} elseif ($background == "shuffle") {
	$background = " -shuffle $shuffle_n,$shuffle_m ";
} elseif ($background == "bgfile") {
	$background = " -bgfile pka.background.txt ";
} elseif ($background == "markov_background") {
	$background = " -markov $markov_background_order -bgfile pka.background.txt ";
} 

	
$startPos = $_POST['startPos'];
$colorblind = $_POST['colorblind'];



// clean up folders older than 3 days
$clean = "find ./files_to_be_removed_72_hrs_after_creation -mtime +3 -exec rm -rf {} \;";
exec('nohup '. $clean . ' > /dev/null 2>&1 &');


// make a folder with randome name for each job
$jobID = substr(str_shuffle(md5(time())),0,20);
 
$tmpfolder = "./files_to_be_removed_72_hrs_after_creation/".$jobID;

$oldmask = umask(0);
if (!mkdir($tmpfolder, 0777, true)) {
    die('Failed to create folders...');
}

// enter this folder
chdir($tmpfolder);

// input
if (strlen($foregroundpaste) > 0) // sequence pasted, save to file
{
	$h = fopen("pka.input.txt", 'w');
	fwrite($h, $foregroundpaste);
	fclose($h); 
}
else
{
	$file1 = $_FILES['foregroundfile'];
	if(move_uploaded_file($file1['tmp_name'], "pka.input.txt")) 
	{
		//echo "The file ".basename( $_FILES['uploadedfile']['name'])." has been uploaded";
	}
	else
	{
		echo "<font color='red'>Please paste or upload your input!</font>";
		exit();
	}
}

// background
if ($background == "bgfile" || $background == "markov_background") {
	if (strlen($backgroundpaste) > 0) // sequence pasted, save to file
	{
		$hb = fopen("pka.background.txt", 'w');
		fwrite($hb, $backgroundpaste);
		fclose($hb); 
	}
	else 
	{
		$file2 = $FILES['backgroundfile'];
		if(move_uploaded_file($file2['tmp_name'], "pka.background.txt"))
		{
		}
		else 
		{
			//echo "<font color='red'>Please paste or upload your background!</font>";
			//exit();
		}
	}
}

//

$command = "PKA pka.input.txt -o pka.output $inputtype -seq $col_seq -weight $col_weight -alphabet $alphabet  $kmer_length -shift $shift $background -startPos $startPos $degenerate $colorblind";

$result = exec('nohup ../../'. $command . ' >> log 2>&1 &');
//$result = exec('nohup ../../'. $command . ' > /dev/null 2>&1 &');

//$myFile = "output.html"; // or .php   
//$fh = fopen($myFile, 'w'); // or die("error");  
//$stringData = "you html code php code goes here";   
//fwrite($fh, $stringData);
//fclose($fh);

umask($oldmask);

exec("cp ../../summary.php ./");

$_SESSION['email'] = $email;
$_SESSION['jobname'] = $jobname;
$_SESSION['command'] = $command;
$_SESSION[$jobID] = "submitted"; // to keep track of job status: submitted, submission_email_sent, results_available_email_sent, 
header("Location: $tmpfolder/summary.php");




?>



