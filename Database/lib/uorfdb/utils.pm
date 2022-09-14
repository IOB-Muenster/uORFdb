# AUTHOR

# Norbert Grundmann


# COPYRIGHT
# THE SOFTWARE IS PROVIDED "AS IS" AND THE AUTHOR DISCLAIMS ALL WARRANTIES WITH
# REGARD TO THIS SOFTWARE INCLUDING ALL IMPLIED WARRANTIES OF MERCHANTABILITY AND
# FITNESS. IN NO EVENT SHALL THE AUTHOR BE LIABLE FOR ANY SPECIAL, DIRECT,
# INDIRECT, OR CONSEQUENTIAL DAMAGES OR ANY DAMAGES WHATSOEVER RESULTING FROM
# LOSS OF USE, DATA OR PROFITS, WHETHER IN AN ACTION OF CONTRACT, NEGLIGENCE
# OR OTHER TORTIOUS ACTION, ARISING OUT OF OR IN CONNECTION WITH THE USE OR
# PERFORMANCE OF THIS SOFTWARE.


package uorfdb::utils;


#
# ===============================================================================
# Basic function for database handling
# ===============================================================================
#

#
# ------------------------------------------------------------------------
# Convert displayed/entered value/state to sql value
# ------------------------------------------------------------------------
#
sub toSQL {
	my ($value, $type) = @_;

	return "null" if (!defined $value);

	if ($type eq "ip") {
		my @a = split(/\./, $value);
		return 0 if (@a != 4);
		return int((($a[3] * 256 + $a[2]) * 256 + $a[1]) * 256 + $a[0]);
	}
	return "$value";
}

1;
