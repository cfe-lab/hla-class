<?php
	session_start();

	if ($_SESSION["results"]) {
		unlink($_SESSION["results"]);
		unlink($_SESSION["details"]);
		unlink($_SESSION["tech_errors"]);
		unlink($_SESSION["errors"]);
		posix_kill($_SESSION["pid"], 15); // SIGTERM, I hope.
	}

	if ($_SESSION["updates"]) {
		unlink($_SESSION["updates"]);
	}

	$ruby = "/usr/local/bin/ruby";
	$script = "hla-easy.rb";
	$_SESSION["results"] = tempnam("./tmp", "");
	$_SESSION["details"] = tempnam("./tmp", "");
	$_SESSION["tech_errors"] = tempnam("./tmp", "");
	$_SESSION["errors"] = tempnam("./tmp", "");

	$letter = $_POST["letter"];
	$descriptors = array (
		0 => array("pipe", "r"),
		1 => array("file", $_SESSION["results"], "a"),
		2 => array("file", $_SESSION["tech_errors"], "a"),
		3 => array("file", $_SESSION["errors"], "a"),
		4 => array("file", $_SESSION["details"], "a")
	);

	$threshold = $_POST["threshold"];

	if ($threshold >= 0) {
		$cmd = sprintf("echo \"%s\" | %s %s -t %s %s &", $_POST["fasta_text"], $ruby,
			$script, $threshold, $letter); 
	} else {
		$cmd = sprintf("echo \"%s\" | %s %s --threshold=\"-1\" %s &", $_POST["fasta_text"], $ruby,
			$script, $letter); 
	}
	$process = proc_open($cmd, $descriptors, $pipes);
	$pstatus = proc_get_status($process);
	$_SESSION["pid"] = $pstatus["pid"];
	fclose($pipes[0]);
	proc_close($process);

	$resfile = "tmp/" . array_pop(explode("/", $_SESSION["results"]));
	$detfile = "tmp/" . array_pop(explode("/", $_SESSION["details"]));
	$techerrfile = "tmp/" . array_pop(explode("/", $_SESSION["tech_errors"]));
	$errfile = "tmp/" . array_pop(explode("/", $_SESSION["errors"]));
	printf("%s\n%s\n%s\n%s", $resfile, $detfile, $techerrfile, $errfile);
?>
