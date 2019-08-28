<!DOCTYPE html>
<html>
	<?php 
		ini_set("error_reporting","E_ALL & ~E_NOTICE"); 
		session_start();
	?>

	<head>
		<title>HLA Classification - New</title>
		<link
		href='https://fonts.googleapis.com/css?family=Numans|Ubuntu+Mono:400,700' rel='stylesheet' type='text/css'>
		<link rel="stylesheet" type="text/css" 
			href="style.css?id=<?php echo rand(1,1000); ?>" />
		<meta name="robots" content="noindex" />  <!-- No thanks search engines! -->
		<script src="scripts.js" type="text/javascript"></script>
	</head>

	<body>
		<div id="tabs">
			<span id="batch_tab" onclick="switchMode('batch')">
				Batch Mode
			</span>
			<span id="single_tab" onclick="switchMode('single')" class="tab_selected">
				Single Sequence Mode
			</span>
			<span id="update_tab" onclick="switchMode('update')">
				Update Allele List
			</span>
		</div>

		<div id="main">
			<!-- Stuff common to both modes. -->
			<form id="main_form" enctype="multipart/form-data" method="post">
			<div id="input">
				<h2>Sequence Input</h2>
				<label>HLA Locus:</label>
				<select id="letter" onchange="toggleOptions(this)">
					<option value="A">A</option>
					<option value="B">B</option>
					<option value="C">C</option>
				</select>
			</div>
	
			<!-- Batch mode only. -->
			<span id="batch_mode" class="hidden">
				<br/><br/>
				<hr>
				<p>Select your files or paste sequences in FASTA format.</p>
		
				<input type="file" id="fasta_file" onchange="displayFile(this.files, 'fasta_text')" />
				<br /><br />
				<textarea id="fasta_text" name="fasta_text" placeholder="Paste your sequences here."></textarea>
			</span>

			<!-- Single mode only. -->
			<span id="single_mode">

				<span id="input_type_choice" class="hidden">
					<label>Sequence type:</label>
					<select id="input_type" onchange="toggleInputType(this.value)">
						<option value="single">Single sequence</option>
						<option value="linked">Exons 2 and 3 separately</option>
					</select>
				</span>
				<br/>
				<label>Maximum number of mismatches:</label>
				<input type="text" id="threshold" value="0" onchange="validateMismatches(this.value)"/>

				<br/><br/>
				<hr>

				<div id="single_mode_single">
					<p>Select sequence file or paste sequence.</p>
					<input type="file" id="fullseq_file" onchange="displayFile(this.files, 'fullseq')" />
					<input type="checkbox" id="fullseq_short" value="true" />Short
					<br /><br />
					<textarea id="fullseq" name="fullseq" placeholder="Paste full sequence here."></textarea>
				</div>

				<div id="single_mode_linked" class="hidden">
					<p>Select exon 2 file or paste sequence.</p>
					<input type="file" id="exon2_file" onchange="displayFile(this.files, 'exon2')" />
					<input type="checkbox" id="exon2_short" value="true" />Short
					<br /><br />
					<textarea id="exon2" name="exon2" placeholder="Paste exon 2 here."></textarea>
					<br /><br />
					<hr>
					<p>Select exon 3 file or paste sequence.</p>
					<input type="file" id="exon3_file" onchange="displayFile(this.files, 'exon3')" multiple/>
					<input type="checkbox" id="exon3_short" value="true" />Short
					<br /><br />
					<textarea id="exon3" name="exon3" placeholder="Paste exon 3 here."></textarea>
				</div>

			</span>

			<span id="update_mode" class="hidden">
				<h2>Update</h2>
				Press "Submit" to update the list of alleles. This will take some time.
				Please do <b>not</b> close the window while the update is in progress,
				otherwise the program will not work correctly.
			</span>

			<br />
			<hr>
			<br />

			<span id="form_errors"></span>
			<button type="button" id="submit" onclick="sendForm()">Submit</button>
			<button type="button" id="reset" onclick="resetForm()">Reset</button>

		</form>

		<span id="status"></span>

		<br/>
		<a id="download" class="hidden" target="_blank">Download csv</a>

		<h2 id="tech_error_header" class="hidden">Technical Errors</h2>
		<div id="tech_error"></div>

		<h2 id="error_header" class="hidden">Errors</h2>
		<div id="error"></div>

		<h2 id="result_header" class="hidden">Results</h2>
		<div id="result"></div>

		<h2 id="detail_header" class="hidden">Possible Matches</h2>
		<div id="detail"></div>

		<h2 id="update_header" class="hidden">Progress</h2>
		<div id="update"></div>
	</body>
</html>
