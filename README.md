## tcolormask. ##

Yet another simple AviSynth plugin created for the sake of creating a simple AviSynth plugin.

This is a port of the great tophf's script [t_colormask][old_script]. It masks colors but unlike the old ColorKeyMask it works on YV12. The main reason for this plugin is performance. The API and output are not identical, please reffer to the **Difference** section for details.

Provided binary is built with ICL. 

### Usage ###
All parameters (except the colors) are set to their default values.
```
tcolormask("$FFFFFF /*pure white hex*/ 
000000 //also black
$808080", tolerance=10, bt601=false, gray=false, lutthr=9, mt=true)
```
* *colors* - string of colors. Line and block comments are allowed.
* *tolerance* - pixel value will pass if its absolute difference with color is less than tolerance (luma) or half the tolerance (chroma).
* *bt601* - use bt601 matrix for conversion of colors.
* *gray* - set chroma of output clip to 128. 
* *lutthr* - if specified more than lutthr colors, lut will be used instead of direct SIMD computations. 
* *mt* - enable multithreading. It actually uses only two threads as any more threading doesn't seem to be useful.

### Difference ###
Unlike the old script, this plugin uses a single string to specify all colors and doesn't do any blurring. Wrapper functions are provided for convenience in the **bin/tcolormask_wrappers.avs** script.

Also since we process chroma and luma together to avoid having to merge planes together later, output is not [identical][comparison]. Basically it appears a bit less blurred and doesn't contain any non-binary values produced by chroma resizing. This will not be fixed and there's no workaround. You'll most likely be blurring the output clip anyway. 

### Performance ###
This plugin uses direct SIMD computations to do its dirty work. SIMD appears to be faster than LUT for the most common cases, but unfortunately its speed depends on the number of specified colors and at some point it does get slower than LUT. That's why the alternative is also provided. If specified more than *lutthr* (9 by default) colors, the plugin will use the LUT routine to avoid performance degradation.

All tests used YV12 1080p image cached by the *loop* function. 

<table>
	<tr>
		<th>
			Colors
		</th>
		<th>
			Calls
		</th>
		<th>
			Old, fps
		</th>
		<th>
			SIMD ST, fps
		</th>
		<th>
			LUT ST, fps
		</th>
		<th>
			SIMD MT, fps
		</th>
		<th>
			LUT MT, fps
		</th>
	</tr>
	<tr>
		<td>
			1
		</td>
		<td>
			1
		</td>
		<td>
			52
		</td>
		<td>
			585
		</td>
		<td>
			150
		</td>
		<td>
			835
		</td>
		<td>
			265
		</td>
	</tr>
	<tr>
		<td>
			1
		</td>
		<td>
			10
		</td>
		<td>
			2
		</td>
		<td>
			82
		</td>
		<td>
			16
		</td>
		<td>
			122
		</td>
		<td>
			31
		</td>
	</tr>
	<tr>
		<td>
			6
		</td>
		<td>
			1
		</td>
		<td>
			52
		</td>
		<td>
			229
		</td>
		<td>
			150
		</td>
		<td>
			396
		</td>
		<td>
			265
		</td>
	</tr>
	<tr>
		<td>
			6
		</td>
		<td>
			10
		</td>
		<td>
			2
		</td>
		<td>
			26
		</td>
		<td>
			16
		</td>
		<td>
			47
		</td>
		<td>
			31
		</td>
	</tr>
	<tr>
		<td>
			24
		</td>
		<td>
			1
		</td>
		<td>
			52
		</td>
		<td>
			71
		</td>
		<td>
			150
		</td>
		<td>
			130
		</td>
		<td>
			265
		</td>
	</tr>
	<tr>
		<td>
			24
		</td>
		<td>
			10
		</td>
		<td>
			2
		</td>
		<td>
			7
		</td>
		<td>
			16
		</td>
		<td>
			14
		</td>
		<td>
			31
		</td>
	</tr>
</table>

### License ###
This project is licensed under the [MIT license][mit_license]. Binaries are [GPL v2][gpl_v2] because if I understand licensing stuff right (please tell me if I don't) they must be.

[old_script]: http://pastebin.com/vAa8fyjp
[comparison]: http://screenshotcomparison.com/comparison/27079
[mit_license]: http://opensource.org/licenses/MIT
[gpl_v2]: http://www.gnu.org/licenses/gpl-2.0.html