<?php
header("Content-Type: octet/stream");
header("Content-Disposition: attachment; filename=".$_GET['name']);
readfile($_GET['path']);
?>