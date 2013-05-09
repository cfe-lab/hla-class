<?php
	session_start();

	posix_kill($_SESSION["pid"], 9);

	unlink($_SESSION["results"]);
	unlink($_SESSION["tech_errors"]);
	unlink($_SESSION["errors"]);
	unlink($_SESSION["details"]);
	unlink($_SESSION["updates"]);

	session_unset();
?>
