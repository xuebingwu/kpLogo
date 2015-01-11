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


// visitor information

$ip = get_real_ip();

$date = date("y.m.d.h.i.s");


$h = fopen("./visitor.txt", 'a');

fwrite($h, $date."\t".$ip."\n");

fclose($h); 

// input	

$inputDir = "./files_to_be_emptied_after_72_hrs/";

$job = $date."-".$ip;

$inputPath = $inputDir.$job;

if (strlen($foreground) > 0) // sequence pasted, save to file
{
	$h = fopen($inputPath, 'w');
	fwrite($h, $foreground);
	fclose($h); 
}
else
{
	if(move_uploaded_file($_FILES['forgroundfile']['tmp_name'], $inputPath)) 
	{
		echo "The file ".basename( $_FILES['uploadedfile']['name'])." has been uploaded";
	}
	else
	{
		echo "<font color='red'>Please paste or upload your input!</font>";
		exit();
	}
}

$outputDir = "./files_to_be_emptied_after_72_hrs/";

$outputPath = $outputDir.$job;

$k = "-k";
if ($kmer_upto == "on")
{
	$k = "-upto";
}

$command = "PKA $inputPath -o $inputPath -alphabet $alphabet $k $kmer_length ";

echo "commandline:"."<br>";

echo $command."<br>";

echo "Input: <a href=\"$inputPath\">download</a> <br>";

echo "Output:<br>";
echo "\t<a href=\"$inputPath\">plot</a> <br>";
echo "\t<a href=\"$inputPath\">data</a> <br>";
echo "\t<a href=\"$inputPath\">summary</a> <br>";



$result = exec($command);

echo $result."<br>";


function get_real_ip()

{

	$ip=false;

	if(!empty($_SERVER["HTTP_CLIENT_IP"]))

	{

		$ip = $_SERVER["HTTP_CLIENT_IP"];

	}

	if (!empty($_SERVER['HTTP_X_FORWARDED_FOR'])) 

	{

		$ips = explode (", ", $_SERVER['HTTP_X_FORWARDED_FOR']);

		if ($ip) 

		{

			array_unshift($ips, $ip); 

			$ip = FALSE;

		}

		for ($i=0;$i<count($ips);$i++)

		{

			if (!eregi ("^(10|172\.16|192\.168)\.", $ips[$i]))

			{

				$ip = $ips[$i];

				break;

			}

		}

	}

	return ($ip ? $ip : $_SERVER['REMOTE_ADDR']);

}



?>
