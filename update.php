<?php
	session_start();

	$ruby = "/usr/local/bin/ruby";
	$script0 = "01_hla_fasta_parse.rb";
	$script1 = "01.5_hla_reduce_by_g.rb";

	$_SESSION["updates"] = tempnam("./tmp", "");
	$_SESSION["tech_errors"] = tempnam("./tmp", "");

	$descriptors = array (
		0 => array("pipe", "r"),
		1 => array("file", $_SESSION["updates"], "w"),
		2 => array("file", $_SESSION["tech_errors"], "w")
	);
	$cmd = sprintf("%s %s && %s %s &", $ruby, $script0, $ruby, $script1);
	$process = proc_open($cmd, $descriptors, $pipes);
	$pstatus = proc_get_status($process);
	$_SESSION["pid"] = $pstatus["pid"];
	fclose($pipes[0]);
	proc_close($process);

	$upfile = "tmp/" . array_pop(explode("/", $_SESSION["updates"]));
	$errfile = "tmp/" . array_pop(explode("/", $_SESSION["tech_errors"]));
	printf("%s\n%s", $upfile, $errfile);
?>
