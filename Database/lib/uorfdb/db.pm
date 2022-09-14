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

package uorfdb::db;

#
# ===============================================================================
# Basic function for database handling
# ===============================================================================
#

use DBI;
use IO::Handle;
use Sys::Hostname;

use lib '../';
use uorfdb::config;

#
# --------------------------------------------------------------------------------
# Returns a handle to the database
# --------------------------------------------------------------------------------
#
sub connect {
	my $dbh = DBI->connect("dbi:Pg:dbname=$DB;host=$HOST", $USER, $PASSWD) or die "Cannot connect to uorf database";
	return $dbh;
}

1;
