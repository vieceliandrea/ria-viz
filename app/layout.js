import './globals.css';

export const metadata = {
  title: 'RIA transient: primary circuit',
  description: 'Reactivity-initiated accident seen from the primary circuit: point kinetics + 1-D thermal hydraulics.',
};
export const viewport = { width: 'device-width', initialScale: 1, viewportFit: 'cover' };

export default function RootLayout({ children }) {
  return (
    <html lang="en" suppressHydrationWarning>
      <head>
        <link rel="preconnect" href="https://fonts.googleapis.com" />
        <link rel="preconnect" href="https://fonts.gstatic.com" crossOrigin="" />
        <link href="https://fonts.googleapis.com/css2?family=Bricolage+Grotesque:wght@600;800&family=Hanken+Grotesk:wght@400;600&display=swap" rel="stylesheet" />
      </head>
      <body>{children}</body>
    </html>
  );
}
