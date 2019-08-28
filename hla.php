<?php
  ini_set("error_reporting","E_ALL & ~E_NOTICE"); 
	session_start();

  //hmmm  
  if ($_SESSION["results"]) {
    unlink($_SESSION["results"]);
    unlink($_SESSION["details"]);
    unlink($_SESSION["tech_errors"]);
    unlink($_SESSION["errors"]);
    //posix_kill($_SESSION["pid"], 15); // SIGTERM, I hope.
	}
    
	if ($_SESSION["updates"]) {
		unlink($_SESSION["updates"]);
	}

	$ruby = "ruby";
	$script = dirname(__FILE__) . "/hla-easy.rb";
	$_SESSION["results"] = dirname(__FILE__) . tempnam("/tmp", "");
	$_SESSION["details"] = dirname(__FILE__) . tempnam("/tmp", "");
	$_SESSION["tech_errors"] = dirname(__FILE__) . tempnam("/tmp", "");
	$_SESSION["errors"] = dirname(__FILE__) . tempnam("/tmp", "");

	$letter = $_POST["letter"];
  $threshold = $_POST["threshold"];
  $fasta_text = $_POST["fasta_text"];
  
  
	$descriptors = array (
		0 => array("pipe", "r"),
		1 => array("file", $_SESSION["results"], "a"),
		2 => array("file", $_SESSION["tech_errors"], "a"),
		7 => array("file", $_SESSION["errors"], "a"),
		8 => array("file", $_SESSION["details"], "a")
	);

	
  
	if ($threshold >= 0) {
		$cmd = sprintf("echo \"%s\" | %s %s -t %s %s &", $fasta_text, $ruby,
			$script, $threshold, $letter); 
	} else {
		$cmd = sprintf("echo \"%s\" | %s %s --threshold=\"-1\" %s &", $fasta_text, $ruby,
			$script, $letter); 
	}
  
	$process = proc_open($cmd, $descriptors, $pipes);
  $pstatus = proc_get_status($process);
  
	$_SESSION["pid"] = $pstatus["pid"];
	fclose($pipes[0]);
	proc_close($process);
  
	$resfile =  "/django/tools/hla_class/tmp/" . array_pop(explode("/", $_SESSION["results"]));
	$detfile =  "/django/tools/hla_class/tmp/" . array_pop(explode("/", $_SESSION["details"]));
	$techerrfile = "/django/tools/hla_class/tmp/" . array_pop(explode("/", $_SESSION["tech_errors"]));
	$errfile = "/django/tools/hla_class/tmp/" . array_pop(explode("/", $_SESSION["errors"]));
	printf("%s\n%s\n%s\n%s", $resfile, $detfile, $techerrfile, $errfile);
?>
