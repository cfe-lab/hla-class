<?php
	session_start();
	
	$lines = file($_SESSION["results"]);
	array_pop($lines);
	file_put_contents($_SESSION["results"], implode("", $lines));
?>
