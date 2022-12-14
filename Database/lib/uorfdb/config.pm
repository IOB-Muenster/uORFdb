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


package uorfdb::config;

use Exporter;
use base 'Exporter';

# Database parameters
our $DB     = "FOO";
our $USER   = "BAR";
our $PASSWD = "FOOBAR";
our $HOST   = "BARFOO";

our @EXPORT = qw($DB $USER $PASSWD $HOST);

1;
