// https://stackoverflow.com/questions/12504042/what-is-a-method-that-can-be-used-to-increment-letters/34483399
function nextChar(c) {
	// increment character.  after 'z' comes 'za', and after 'za' comes 'zb'
	
	// get last character of string to increment
	if (c.length > 1) {
		var toIncrement = c.substring(c.length-1, c.length);
		var restOfString = c.substring(0, c.length-1);
	} else {
		var toIncrement = c;
		var restOfString = '';
	}
	
	// if last character is z or Z, add a to end of string and return
	if (toIncrement == 'z') {
		return 	c.concat('a');
	} else if (toIncrement == 'Z') {
		return c.concat('A');
	} else {
		return restOfString.concat(String.fromCharCode(toIncrement.charCodeAt(0) + 1));	
	}
}

// https://www.javascripttutorial.net/dom/manipulating/insert-an-element-after-an-existing-element/
function insertAfter(newNode, existingNode) {
	existingNode.parentNode.insertBefore(newNode, existingNode.nextSibling);
}	


// create a new set of barcodes
function createNewSet(type) {
	
	// get all group fieldsets
	const fieldsets = document.getElementsByClassName("group set");

	// get letter for this fieldset	
	if (fieldsets.length == 0) {
		var next_let = "A";
	} else {
		var last_set = fieldsets[fieldsets.length - 1].getAttribute("id").split("_")[0];
		var next_let = nextChar(last_set);
	}

	console.log(`Creating new ${type} set ${next_let}`);
	
	// get set number
	const set_num = next_let.charCodeAt(0) - 64;
	
	// create new fieldset
	var new_set = document.createElement("fieldset");
	new_set.setAttribute("id", `${next_let}_set`);
	new_set.setAttribute("class", `${next_let} set group`);

	// create a label for the set
	var label = document.createElement("p");
	label.setAttribute("class", "`${set_let} set`");
	var text = document.createTextNode(`Set ${set_num} (${type}) name: `);
	label.appendChild(text);
	
	// create a name box for the set
	var input = document.createElement("input");
	input.setAttribute("type", "text");
	input.setAttribute("id", `set_${next_let}_name`);
	input.setAttribute("class", `${next_let} set name`);
	input.setAttribute("name", `set_${next_let}_name`);
	input.setAttribute("placeholder", "Name");
	
	// append name box to label paragraph
	label.append(input)
	
	// append label as a child of new_set
	new_set.appendChild(label);
	
	return {
		next_let,
		new_set
	};

}

// add new barcode set at end of other sets
function addSet(new_set) {

	var buttons = document.getElementsByTagName("button");
	var second_last_button = buttons[buttons.length-2];
	
	second_last_button.parentNode.insertBefore(new_set, second_last_button);

}


// Add a new constant set
function addConstSet() {
	
	// create new fieldset
	var new_set = createNewSet("constant");
	
	// add position input box
	var label = document.createElement("p");
	label.setAttribute("class", "`${new_set.next_let} set pos`");
	var text = document.createTextNode(`Position in read: `);
	label.appendChild(text);
	
	// create a name box for the set
	var input = document.createElement("input");
	input.setAttribute("type", "text");
	input.setAttribute("id", `set_${new_set.next_let}_pos`);
	input.setAttribute("class", `${new_set.next_let} pos`);
	input.setAttribute("name", `set_${new_set.next_let}_pos`);
	input.setAttribute("placeholder", "0-based position");
	
	label.append(input);
	new_set.new_set.appendChild(label);
	
	// add the first barcode
	var barc = createBarcode(new_set.next_let, "a");
	new_set.new_set.append(barc);
	
	// add a button for adding more barcodes
	var but = document.createElement("button");
	but.setAttribute("id", `${new_set.next_let}_more`);
	but.setAttribute("onClick", "addBarcode(this.id)");
	but.setAttribute("type", "button");
	var but_text = document.createTextNode("Add barcode");
	but.appendChild(but_text);
	
	new_set.new_set.appendChild(but);
	
	addSet(new_set.new_set);
	
}

// Add a new variable set of barcodes
function addVarSet() {

	// create new fieldset
	var new_set = createNewSet("variable");
	
	// add a 'before' input box
	var label = document.createElement("p");
	label.setAttribute("class", "`${new_set.next_let} set before`");
	var text = document.createTextNode(`Before sequence: `);
	label.appendChild(text);
	
	// create a name box for the set
	var input = document.createElement("input");
	input.setAttribute("type", "text");
	input.setAttribute("id", `set_${new_set.next_let}_before`);
	input.setAttribute("class", `${new_set.next_let} before`);
	input.setAttribute("name", `set_${new_set.next_let}_before`);
	input.setAttribute("placeholder", "Sequence");
	
	label.append(input);
	new_set.new_set.appendChild(label);	
	
	// add an 'after' input box
	var label = document.createElement("p");
	label.setAttribute("class", "`${new_set.next_let} set after`");
	var text = document.createTextNode(`After sequence: `);
	label.appendChild(text);
	
	// create a name box for the set
	var input = document.createElement("input");
	input.setAttribute("type", "text");
	input.setAttribute("id", `set_${new_set.next_let}_after`);
	input.setAttribute("class", `${new_set.next_let} after`);
	input.setAttribute("name", `set_${new_set.next_let}_after`);
	input.setAttribute("placeholder", "Sequence");
	
	label.append(input);
	new_set.new_set.appendChild(label);	
	
	// add a check box to say if we should translate these variable barcodes or not
	var box_label = document.createElement("label");
	box_label.setAttribute("for", `${new_set.next_let}_translate`);
	var box_text = document.createTextNode(`Translate extracted sequences: `);
	box_label.appendChild(box_text);
	new_set.new_set.append(box_label)

	var box = document.createElement("input");
	box.setAttribute("type", "checkbox");
	box.setAttribute("id", `${new_set.next_let}_translate`);
	box.setAttribute("name", `${new_set.next_let}_translate`);
	box.setAttribute("value", "Translate");
	new_set.new_set.appendChild(box)

	addSet(new_set.new_set);

}
// function to create a new barcode fieldset
function createBarcode(set_let, next_let) {

	// create a new fieldset
	var new_set = document.createElement("fieldset");
	new_set.setAttribute("id", [set_let, next_let, "barc"].join("_"));
	new_set.setAttribute("class", "`${set_let} ${next_let} barc field`")
	
	// create a label for the barcode
	var label = document.createElement("p");
	label.setAttribute("class", "`${set_let} ${next_let} barc`")
	var text = document.createTextNode("Barcode: ");
	label.appendChild(text);
	
	// create two input elements
	var bname = document.createElement("input");
	bname.setAttribute("type", "text");
	bname.setAttribute("name", [set_let, next_let, "name"].join("_"));
	bname.setAttribute("placeholder", "Name");
	bname.setAttribute("class", `${set_let} ${next_let} name barc`);
	bname.setAttribute("id", `${set_let}_${next_let}_name`);
	
	var bseq = document.createElement("input");
	bseq.setAttribute("type", "text");
	bseq.setAttribute("name", [set_let, next_let, "seq"].join("_"));
	bseq.setAttribute("placeholder", "Sequence");
	bseq.setAttribute("class", `${set_let} ${next_let} seq barc`);
	bseq.setAttribute("id", `${set_let}_${next_let}_seq`);
	
	label.appendChild(bname);
	label.appendChild(bseq);
	
	new_set.append(label);
	
	return new_set;		
}

// function to add a new barcode fieldset when button with id is clicked
function addBarcode(id) {

	console.log(id)

	// get set letter            
	var info = id.split("_");
	var set_let = info[0];
	
	console.log(set_let)
	
	// get letter for last barcode in the set
	var inputs = document.querySelectorAll(`input.${set_let}`);
	
	var last_let = inputs[inputs.length-1].getAttribute("id").split("_")[1];
	
	// increment for name of next barcode
	var next_let = nextChar(last_let);
	
	console.log(`Adding new barcode ${next_let} to set ${set_let}`);
	
	new_set = createBarcode(set_let, next_let);
	
	console.log(new_set);
	
	insertAfter(new_set, document.querySelector(`#${set_let}_${last_let}_barc`));
	
}


function validateForm() {

	return true;
}

