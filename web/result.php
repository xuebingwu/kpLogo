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
		echo "The file ".basename( $_FILES['uploadedfile']['name'])." has been uploaded";
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

$command = "/lab/bartel1_ata/wuxbl/scripts/C++/bin/PKA PKA.input.txt -o PKA";

echo "commandline:"."<br>";

echo $command."<br>";


$result = exec($command);

umask($oldmask);

echo "Input: <a href=\"$tmpfolder/PKA.input.txt\">download</a> <br>";

echo "Output (to be removed after 72 hours):<br>";
echo "\t<a href=\"$tmpfolder/PKA.freq.ps\">frequency logo</a> <br>";
echo "\t<a href=\"$tmpfolder/PKA.info.ps\">information content logo</a> <br>";
echo "\t<a href=\"$tmpfolder/PKA.ps\">significance logo</a> <br>";
echo "\t<a href=\"$tmpfolder/PKA.most.significant.each.position.ps\">kmer logo</a> <br>";
echo "\t<a href=\"$tmpfolder/PKA.most.significant.each.position.txt\">most.significant.kmer.at.each.position</a> <br>";
echo "\t<a href=\"$tmpfolder/PKA.pass.p.cutoff.txt\">output data file</a> <br>";

echo "done <br>";


?>
